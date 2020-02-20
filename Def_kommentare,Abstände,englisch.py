#!/usr/bin/env python3
import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
import time
#from alive_progress import alive_bar #Process bar
#from alve_progress import config_handler #Process bar
from itertools import zip_longest, chain
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec


class Read:
    """Store read with all necessary information as attributes."""
    def __init__(self, content, id, phred, alphabet, orig_quality, trim_quality):
        self.length = len(content[1])
        self.gc = (content[1].count("G") + content[1].count("C"))/self.length #in decimal - between 0 and 1
        self.name = content[0].split()[0] #ID taken from firt line
        self.orig_quality = orig_quality
        self.quality = trim_quality
        self.sequenz = content[1][0:len(self.quality)] #sequence-string
        self.mean_qual_trimmed = np.mean(trim_quality)
        self.mean_qual_untrimmed = np.mean(orig_quality)
        self.trimmed_length = len(trim_quality)


def parser():
    """Parse command line input when starting script."""
    input_parser = argparse.ArgumentParser(description="Quality evaluation for FASTA files")
    input_parser.add_argument(
        "File path",
        help="File path to the FASTA file"
        )
    input_parser.add_argument(
        "-p", "--phred",
        dest="phred",
        default=None,
        help="Optional manual specification of the Phred format, 33 or 64"
        )
    input_parser.add_argument(
        "-t", "--trim",
        dest="trim_val",
        help="Optional manual specification of the trimming score",
        default=25
        )
    input_parser.add_argument(
        "-c", "--cutoff",
        dest="cutoff",
        help="Optional manual specification of the average score,where a read is completely discarded.",
        default=20
        )
    input_parser.add_argument(
        "-m", "--minl",
        dest="minlength",
        help="Optional manual specification of the minimum length in bases that a trimmed read must have during evaluation.",
        default=20
        )
    input_parser.add_argument(
        "-tab", "--table",
        dest="save_table",
        help="Optional specification of a file path for saving the data table.",
        default="Datentabelle.csv"
        )
    input_parser.add_argument(
        "-plt", "--plot",
        dest="save_plot",
        help="Optional specification of a file path for saving the data plots.",
        )
    input_parser.add_argument(
        "-i", "--interaktiv",
        dest="interaktiv",
        action="store_true",
        help="Activation of the commented interactive mode.",
        )
    arguments = input_parser.parse_args() #Creates a dictionary, access to argument by name

    return arguments


def alphabet_ascii(phred, ascii):
    """Create dictionary with scores for ASCII by given phred format."""
    alphabet = None

    if phred == "33": #check for phred format
        alphabet = {char:value-33 for value,char in enumerate(ascii,33) if value < 74} #make dictionary w ascii-character + phred score

    elif phred == "64":
        alphabet = {char:value-64 for value,char in enumerate(ascii,33) if value > 63 and value < 105} #make dictionary w ascii-character + phred score

    return alphabet


def phred_bestimmung(datei, ascii, phred=None):
    """Check for manually given phred.
    If not given, try to get phred from first line.
    If no unique ASCII found, quit program.
    """
    if not phred:

        for score_char in datei:

            if score_char in [char for value, char in enumerate(ascii, 33) if value < 64]:
                phred = "33"
                break

            elif score_char in [char for value, char in enumerate(ascii, 33) if value > 63 and value < 105]:
                phred = "64"
                break

            else:
                print("Phred-Score nicht bestimmbar. Bitte manuell mit -p eingeben.")
                quit()

    return phred


def qualitaet(content, alphabet, phred):
    """Create Phred-Score dictionary and translate ASCII-character into score.
    Append score to value list.
    """
    values = []
    values += [alphabet[str(letter)] for letter in content[3]]

    return values


def trimming(scores, trim_val):
    """Trim sequence and quality list by given trim value.
    Calculate mean score of kmers and cuts off end if lower value found.
    """
    counter = 0
    summe = 0
    item_index = 0

    for kmer in scores:
        counter += 1
        summe += kmer
        item_index += 1

        if counter == 3:
            kmer_mean = summe / 3
            counter = 0
            summe = 0

            if kmer_mean < trim_val:
                return scores[0:item_index]

    else:
        return scores


def graph_gc(dict_reads):
    """Create subplot for gc-percentages per read."""
    fig = plt.figure()
    ax = fig.subplots(1, 1)
    fig.subplots_adjust(hspace=0.5)
    ax.set(
        title="GC Content",
        xlabel="GC Content in %",
        ylabel="Number of Reads",
        )
    sns.distplot(
        [int(round(item.gc,2)*100) for item in dict_reads],
        bins=20,
        kde=False,
        ax=ax,
        color="teal",
        hist_kws=({"alpha":1}),
        )

    return fig


def graph_len_count(dict_reads):
    """Create graph showing read length distribution
    over the reads bevore and after trimming."""
    maxLen = max([item.length for item in dict_reads])

    fig = plt.figure()
    ax1, ax2 = fig.subplots(2, 1, sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.2)
    ax1.set(
        title="Readlengths before Trimming",
        ylabel="Number of Reads",
        xlim=(0, maxLen)
        )
    ax1.set_xlim(0, maxLen)
    sns.distplot(
        [item.length for item in dict_reads],
        bins=40,
        kde=False,
        ax=ax1,
        color="teal",
        hist_kws=({"alpha":1}),
        )

    ax2.set(
    title="Readlengths after Trimming",
        xlabel="Readlength in BP",
        ylabel="Number of Reads",
        )
    
    sns.distplot(
        [item.trimmed_length for item in dict_reads],
        bins=40,
        kde=False,
        ax=ax2,
        color="teal",
        hist_kws=({"alpha":1}),
        )

    return fig


def graph_scores(dict_reads):
    """Create subplot for mean scores per read. Shows the average quality
    of the reads bevore and after trimming """
    fig = plt.figure()
    ax1, ax2 = fig.subplots(2, 1, sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.2)

    ax1.set(
        title="Average Sequence Quality before Trimming",
        ylabel="Number of Reads",
        xlim=(0,40),
        )
    
    sns.distplot(
    [item.mean_qual_untrimmed for item in dict_reads],
        bins=40,
        kde=False,
        ax=ax1,
        color="teal",
        hist_kws=({"alpha":1})
        )

    ax2.set(
        title="Average Sequence Quality after Trimming",
        xlabel="Quality Score",
        ylabel="Number of Reads",
        )
    
    sns.distplot(
    [item.mean_qual_trimmed for item in dict_reads],
        bins=40,
        kde=False,
        ax=ax2,
        color="teal",
        hist_kws=({"alpha":1})
        )

    return fig


def graph_basequality(dict_reads):
    """Show the average quality of all read positions including the standard deviation"""
    qualscores = [item.quality for item in dict_reads]
    df_qual = pd.DataFrame(list(zip_longest(*qualscores,fillvalue=np.nan)), dtype=float)
    df_qual["Sterr"] = [np.nanstd(item) for item in df_qual.values]
    df_qual["Mean"] = [np.nanmean(item) for item in df_qual.values]

    fig = plt.figure()
    ax = fig.subplots(1, 1, sharex=True, sharey=True)
    ax.set(
        title="Average Quality at Read Positions",
        xlabel="Position in Read (BP)",
        ylabel="Quality",
        ylim=(0, df_qual["Mean"].max()+df_qual["Sterr"].max()),
        xlim=(0, df_qual.shape[0]),
        )
    
    plt.errorbar(
        x=df_qual.index.values.tolist(),
        y=df_qual["Mean"].tolist(),
        yerr=df_qual["Sterr"].tolist(),
        data=df_qual,
        color="#008080",
        lw=0.8,
        ecolor="lightblue",
        capthick="0.5",
        )

    return fig


def graph_basenvert(dict_reads):
    """Show average percentage of each bases at read positions"""
    sequences = [list(item.sequenz) for item in dict_reads]
    datalist = list(zip_longest(*sequences, fillvalue=""))
    bases = ["A", "C", "G", "T", "N"]
    colors=["red", "green", "blue", "deeppink", "purple"]
    counter = 0

    fig = plt.figure()
    ax = fig.subplots(1,1)
    
    for letter in bases:
        mean=[i.count(letter)/len("".join(i))*100 for i in datalist]
        y_values=[np.nanmean(pack) for pack in [mean[teil:teil+5] for teil in range(0, len(mean), 5)]]
        sns.lineplot(
                x=[index*5 for index in range(0,max([len(y_values)]))],
                y=y_values,
                lw=0.7,
                color=colors[counter],
                )
        counter += 1

    plt.title(f"Average Percentage of Bases at Read Positions")
    plt.xlabel("Position in Read (BP)")
    plt.ylabel(f"Share of Bases over all Reads (%)")
    plt.xlim(0,max(len(item) for item in sequences))
    plt.ylim(0,100)
    
    fig.legend(
        bases,
        loc="upper center",
        ncol=5,
        fancybox=True,
        framealpha=1,
        edgecolor="#000000",
        facecolor="#FFFFFF",
        )
    
    return fig


def graph_basenvert_abr(dict_reads):
    """Show average percentage of each bases in packs of 5 at read positions
    and in enlarged presentation"""
    fig = plt.figure()
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.5, wspace=0.3)
    gs = GridSpec(2, 2)
    sequences = [list(item.sequenz) for item in dict_reads]
    datalist = list(zip_longest(*sequences, fillvalue=""))
    bases = ["A", "C", "G", "T"]
    colors=["red", "green", "blue", "deeppink"]
    counter = 0

    for letter in bases:
        mean=[i.count(letter)/len("".join(i))*100 for i in datalist]
        y_values=[np.nanmean(pack) for pack in [mean[teil:teil+5] for teil in range(0, len(mean), 5)]]
        ax = fig.add_subplot(list(gs)[counter])
        sns.lineplot(
                x=[index*5 for index in range(0,max([len(y_values)]))],
                y=y_values,
                lw=0.7,
                color=colors[counter],
                )

        plt.title(f"Distribution of {letter}", fontdict={"size":10})
        plt.xlabel("Position in Read (BP)",fontdict={"size":8})
        plt.ylabel(f"Share of {letter} (%)",fontdict={"size":8})
        counter += 1
    plt.suptitle(f"Average Percentage of Bases at Read Positions")

    return fig


def graph_nshare(dict_reads):
    '''Show average percentage of each bases in packs of 5 at read positions
    and in enlarged presentation'''
    sequences = [list(item.sequenz) for item in dict_reads]
    datalist = list(zip_longest(*sequences, fillvalue=""))
    mean=[i.count("N")/len("".join(i))*100 for i in datalist]
    fig = plt.figure()
    ax = fig.subplots(1,1)
    y_values=[np.nanmean(pack) for pack in [mean[teil:teil+3] for teil in range(0, len(mean), 3)]]
    
    sns.lineplot(
            x=[index*5 for index in range(0,max([len(y_values)]))],
            y=y_values,
            lw=0.7,
            color="purple",
            )
    
    plt.title(f"Average Percentage of N at Read Positions")
    plt.xlabel("Position in Read (BP)")
    plt.ylabel(f"Share of N over all Reads (%)")

    return fig


def tabelle_speichern(dict_reads, phred, read_count, dateipfad="Datentabelle.csv"):
    """Save data stats/attributes as csv table for later review. Argument has a default"""
    orig_length, trim_length, orig_mean_qual, trim_mean_qual, gc_con = [[],[],[],[],[]]

    for item in dict_reads:
        orig_length.append(item.length)
        trim_length.append(item.trimmed_length)
        orig_mean_qual.append(item.mean_qual_untrimmed)
        trim_mean_qual.append(item.mean_qual_trimmed)
        gc_con.append(item.gc*100)

    dataframe_overview = pd.DataFrame({
        "Anzahl Reads":[
            read_count,
            "",
            len(dict_reads),
            ""
            ],
        "Phred-Score Sytem":phred,
        "Durchschn. Länge":[
            np.mean(orig_length),
            np.std(orig_length),
            np.mean(trim_length),
            np.std(trim_length),
            ],
        "Durchschn. Qualität":[
            np.mean(orig_mean_qual),
            np.std(orig_mean_qual),
            np.mean(trim_mean_qual),
            np.std(trim_mean_qual),
            ],
        "Durchschn. GC-Gehalt %":[
            "",
            "",
            np.mean(gc_con),
            np.std(gc_con),
            ],
        "":""
        })

    dataframe_overview.index = np.arange(1,len(dataframe_overview)+1)
    
    with open(dateipfad, "w"):
        dataframe_overview.transpose().to_csv(
            path_or_buf=dateipfad,
            header=[
                "Vor Trimming",
                "StAbw vorher",
                "Nach Trimming",
                "StAbw nachher"
                ],
            encoding="utf-8"
            )

    dataframe_detail = pd.DataFrame({
        "Name":item.name,
        "Länge vor Trimming":item.length,
        "Länge nach Trimming":item.trimmed_length,
        "GC":item.gc,
        "Qualität":str(item.quality).strip("[]"),
        "Sequenz":item.sequenz
        } for item in dict_reads)
    dataframe_detail.index = np.arange(1,len(dataframe_detail)+1)

    with open(dateipfad, "a"):
            dataframe_detail.to_csv(
                path_or_buf=dateipfad,
                header=dataframe_detail.columns,
                index_label="Index",
                mode="a",
                encoding="utf-8"
                )


def main():
    #config_handler.set_global(spinner="fish_bouncing", bar="blocks")#Process bar
    #with alive_bar() as bar: #Process bar
    all_ids = [] #! saved as object code, not as string!
    ascii = ("!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
    arguments = parser()
    cutoff = arguments.cutoff
    minlength = arguments.minlength
    phred = arguments.phred
    trim_val = arguments.trim_val

    with open(arguments.Dateipfad) as inhalt:
        block = []
        #test_cnt = 0	#!! KÜRZEN FÜR TESTLÄUFE

        for lines in inhalt:

            if block == []:
                phred = phred_bestimmung(
                    lines,
                    ascii,
                    phred
                    )

            block.append(lines.rstrip())
            #test_cnt += 1	#!! KÜRZEN FÜR TESTLÄUFE

            #if test_cnt > 10000:	#!! KÜRZEN FÜR TESTLÄUFE
            #    break

    alphabet = alphabet_ascii(
        phred,
        ascii
        )
    line_pack = []
    read_count = 0

    for lines in block:
        line_pack.append(lines)

        if len(line_pack) == 4:
            read_count += 1
            quality = qualitaet(
                line_pack,
                alphabet,
                phred
                )
            trim_quality = trimming(
                quality,
                trim_val
                )
            if np.mean(trim_quality) > cutoff and len(trim_quality) >= minlength:
                all_ids.append(
                    Read(
                        line_pack,
                        id,
                        phred,
                        alphabet,
                        quality,
                        trim_quality
                        )
                    )
            line_pack = []
    print("Finished reading file and calculations.")


    #Following code is responsible for creation and saving of graphs
    print("Starting plots.")
    if arguments.save_plot:
        with PdfPages(str(arguments.save_plot)) as pdf:
            #Multipage-pdf with all plots
            sns.set_style("whitegrid", {
                "grid_linestyle":"+",
                "grid.color":"#99ACAC",
                "xtick.bottom":True,
                "ytick.left":True,
                "ytick.direction":"in",
                "axes.edgecolor":"#004A56",
                "axes.facecolor":"F9FEFE",
                })
            
            graph_gc(all_ids)
            pdf.savefig()
            plt.close()

            fig = graph_len_count(all_ids)
            pdf.savefig()
            plt.close()

            graph_scores(all_ids)
            pdf.savefig()
            plt.close()

            graph_basequality(all_ids)
            pdf.savefig()
            plt.close()

            graph_basenvert(all_ids)
            pdf.savefig()
            plt.close()

            graph_basenvert_abr(all_ids)
            pdf.savefig()
            plt.close()

            graph_nshare(all_ids)
            pdf.savefig()
            plt.close()
            end = time.time()
            print(end-start)

    else:
        graph_gc(all_ids)
        plt.show()
        graph_len_count(all_ids)
        plt.show()
        graph_scores(all_ids)
        plt.show()
        graph_basequality(all_ids)
        plt.show()
        graph_basenvert(all_ids)
        plt.show()
        graph_basenvert_abr(all_ids)
        plt.show()
        graph_nshare(all_ids)
        plt.show()

    tabelle_speichern(all_ids, phred, read_count, arguments.save_table)
    print("Table saved.")
    #bar(text="Loading...")#Process bar

if __name__ == "__main__":
    start = time.time()
    main()
    end = time.time()
    print(end-start)
