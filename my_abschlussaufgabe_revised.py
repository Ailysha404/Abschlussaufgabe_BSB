#! /usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
import time
import csv
from alive_progress import alive_bar, config_handler
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
    input_parser = argparse.ArgumentParser(description="Qualitätsauswertung für FASTA-Dateien")
    input_parser.add_argument(
        "Dateipfad",
        help="Dateipfad zur FASTA-Datei"
        )
    input_parser.add_argument(
        "-p", "--phred",
        dest="phred",
        default=None,
        help="Optionale manuelle Angabe des Phred-Formats, 33 oder 64"
        )
    input_parser.add_argument(
        "-t", "--trim",
        dest="trim_val",
        help="Optionale manuelle Angabe des Trimming-Scores",
        default=25
        )
    input_parser.add_argument(
        "-c", "--cutoff",
        dest="cutoff",
        help="Optionale manuelle Angabe des Durchschnittscores, bei dem ein Read komplett verworfen wird.",
        default=20
        )
    input_parser.add_argument(
        "-m", "--minl",
        dest="minlength",
        help="Optionale manuelle Angabe der minimalen Länge in Basen, die ein getrimmter Read bei der Auswertung haben muss.",
        default=20
        )
    input_parser.add_argument(
        "-tab", "--table",
        dest="save_table",
        help="Optionale Angabe eines Dateipfades zum Speichern der Datentabelle",
        default="Datentabelle.csv"
        )
    input_parser.add_argument(
        "-plt", "--plot",
        dest="save_plot",
        help="Optionale Angabe eines Dateipfades zum Speichern der Datenplots",
        )
    arguments = input_parser.parse_args() #Erstellen dictionary, Zugriff auf Argument über Namen

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

def calc_data(dict_reads, phred):
    name, orig_length, trim_length, orig_qual, trim_qual, gc_con, mean_qual_trim_base, sequences = [[],[],[],[],[],[],[],[]]

    qual_pos_sterr = []
    qual_pos_mean = []
    dict_base_pos = {}
    bases = ["A", "C", "G", "T", "N"]
    count = 0
    read_count = len(dict_reads)

    for item in dict_reads:
        name.append(item.name)
        orig_length.append(item.length)
        trim_length.append(item.trimmed_length)
        orig_qual.append(item.mean_qual_untrimmed)
        trim_qual.append(item.mean_qual_trimmed)
        gc_con.append(round(item.gc, 2)*100)
        mean_qual_trim_base.append(item.quality)
        sequences.append(item.sequenz)

    zip_sequenz = list(zip_longest(*sequences, fillvalue=""))
    for letter in bases:
        mean = [i.count(letter)/len("".join(i))*100 for i in zip_sequenz]
        dict_base_pos[letter] = [np.nanmean(pack) for pack in [mean[teil:teil+5] for teil in range(0, len(mean), 5)]]

    for item in list(zip_longest(*mean_qual_trim_base,fillvalue=np.nan)):
        count += 1
        qual_pos_sterr.append(np.nanstd(item))
        qual_pos_mean.append(np.nanmean(item))

    series_data = pd.Series(data={
        "Name":name,
        "Number of Reads":read_count,
        "Phred-Score Sytem":phred,
        "Length before Trimming":orig_length,
        "Avg Length before Trimming":(sum(orig_length)/read_count),
        "Max Length before Trimming":max(orig_length),
        "Length after Trimming":trim_length,
        "Avg Length after Trimming":(sum(trim_length)/read_count),
        "Max Length after Trimming":max(trim_length),
        "Quality before Trimming":orig_qual,
        "Avg Quality before Trimming":(sum(orig_qual)/read_count),
        "Quality after Trimming":trim_qual,
        "Avg Quality after Trimming":(sum(trim_qual)/read_count),
        "GC Share in %": gc_con,
        "Avg GC Share %":(sum(gc_con)/read_count),
        "Mean Quality per Base after Trimming":qual_pos_mean,
        "Max Mean Quality per Base after Trimming":max(qual_pos_mean),
        "Std Error Quality per Base after Trimming":qual_pos_sterr,
        "Max Error Quality per Base after Trimming":max(qual_pos_sterr),
        "Mean A Content per Base after Trimming":dict_base_pos["A"],
        "Mean C Content per Base after Trimming":dict_base_pos["C"],
        "Mean G Content per Base after Trimming":dict_base_pos["G"],
        "Mean T Content per Base after Trimming":dict_base_pos["T"],
        "Mean N Content per Base after Trimming":dict_base_pos["N"],
        "Sequenz":sequences,
        })

    return series_data

def graph_gc(series_data): #opimize
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
        series_data["GC Share in %"],
        bins=20,
        kde=False,
        ax=ax,
        color="teal",
        hist_kws=({"alpha":1}),
        )

    return fig


def graph_len_count(series_data):
    """Create graph showing read length distribution over the reads."""
    maxLen = series_data["Max Length before Trimming"]

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
        series_data["Length before Trimming"],
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
        series_data["Length after Trimming"],
        bins=40,
        kde=False,
        ax=ax2,
        color="teal",
        hist_kws=({"alpha":1}),
        )

    return fig


def graph_scores(series_data):
    """Create subplot for mean scores per read."""
    fig = plt.figure()
    ax1, ax2 = fig.subplots(2, 1, sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.2)

    ax1.set(
        title="Average Sequence Quality before Trimming",
        ylabel="Number of Reads",
        xlim=(0,40),
        )
    sns.distplot(
        series_data["Quality before Trimming"],
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
        series_data["Quality after Trimming"],
        bins=40,
        kde=False,
        ax=ax2,
        color="teal",
        hist_kws=({"alpha":1})
        )

    return fig


def graph_basequality(series_data):
    """Show the quality of all Read positions"""
    fig = plt.figure()
    ax = fig.subplots(1, 1, sharex=True, sharey=True)
    ax.set(
        title="Average Quality at Read Positions",
        xlabel="Position in Read (BP)",
        ylabel="Quality",
        ylim=(0, series_data["Max Mean Quality per Base after Trimming"]+series_data["Max Error Quality per Base after Trimming"]),
        xlim=(0, series_data["Max Length after Trimming"]),
        )
    plt.errorbar(
        x=range(0,int(series_data["Max Length after Trimming"])),
        y=series_data["Mean Quality per Base after Trimming"],
        yerr=series_data["Std Error Quality per Base after Trimming"],
        data=series_data,
        color="#008080",
        lw=0.8,
        ecolor="lightblue",
        capthick="0.5",
        )

    return fig


def graph_basenvert(series_data):
    """Show the distribution of all bases"""
    colors=["red", "green", "blue", "deeppink", "purple"]
    bases = ["A", "C", "G", "T", "N"]
    counter = 0

    fig = plt.figure()
    ax = fig.subplots(1,1)
    for letter in bases:
        sns.lineplot(
                x=range(0,series_data["Max Length after Trimming"], 5),
                y=series_data[f"Mean {letter} Content per Base after Trimming"],
                lw=0.7,
                color=colors[counter],
                )
        counter += 1

    plt.title(f"Average Percentage of Bases at Read Positions")
    plt.xlabel("Position in Read (BP)")
    plt.ylabel(f"Share of Bases over all Reads (%)")
    plt.xlim(0,series_data["Max Length after Trimming"])
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


def graph_basenvert_abr(series_data):
    """Show the distribution of all bases"""
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.5, wspace=0.3)
    gs = GridSpec(2, 2)
    bases = ["A", "C", "G", "T"]
    colors=["red", "green", "blue", "deeppink"]
    counter = 0

    for letter in bases:
        ax = fig.add_subplot(list(gs)[counter])
        sns.lineplot(
                x=range(0,series_data["Max Length after Trimming"], 5),
                y=series_data[f"Mean {letter} Content per Base after Trimming"],
                lw=0.7,
                color=colors[counter],
                )

        plt.title(f"Distribution of {letter}", fontdict={"size":10})
        plt.xlabel("Position in Read (BP)",fontdict={"size":8})
        plt.ylabel(f"Share of {letter} (%)",fontdict={"size":8})
        counter += 1
    plt.suptitle(f"Average Percentage of Bases at Read Positions")

    return fig


def graph_nshare(series_data):
    fig = plt.figure()
    ax = fig.subplots(1,1)
    sns.lineplot(
            x=range(0,series_data["Max Length after Trimming"], 5),
            y=series_data["Mean N Content per Base after Trimming"],
            lw=0.7,
            color="purple",
            )
    plt.title(f"Average Percentage of N at Read Positions")
    plt.xlabel("Position in Read (BP)")
    plt.ylabel(f"Share of N over all Reads (%)")

    return fig

def tabelle_speichern(series_data, dateipfad="Datentabelle.csv"):
    """Save data stats/attributes as csv table for later review"""
    with open(dateipfad, "w"):
        series_data.get(["Number of Reads", "Phred-Score Sytem", "Avg Length before Trimming", "Max Length before Trimming", "Avg Length after Trimming", "Max Length after Trimming", "Avg Quality before Trimming", "Avg Quality after Trimming", "Avg GC Share %"]).to_csv(
            path_or_buf=dateipfad,
            header=False,
            encoding="utf-8",
            )
    detail = series_data.get(["Name", "Quality before Trimming", "Quality after Trimming", "Length before Trimming", "Length after Trimming", "GC Share in %", "Sequenz"])

    mydata=[]
    mydata.append([key for key, item in detail.items()])
    mydata.extend(list(zip(*[item for key, item in detail.items()])))

    with open(dateipfad, "a"):
        for zeile in mydata:
            writer = csv.writer(open(dateipfad, "a"))
            writer.writerow(zeile)


def main():
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

            #if test_cnt > 100:	#!! KÜRZEN FÜR TESTLÄUFE
                #break

    alphabet = alphabet_ascii(
        phred,
        ascii
        )
    line_pack = []

    for lines in block:
        line_pack.append(lines)

        if len(line_pack) == 4:
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

    all_data_series = calc_data(all_ids, phred)
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
            graph_gc(all_data_series)
            pdf.savefig()
            plt.close()

            fig = graph_len_count(all_data_series)
            pdf.savefig()
            plt.close()

            graph_scores(all_data_series)
            pdf.savefig()
            plt.close()

            graph_basequality(all_data_series)
            pdf.savefig()
            plt.close()

            graph_basenvert(all_data_series)
            pdf.savefig()
            plt.close()

            graph_basenvert_abr(all_data_series)
            pdf.savefig()
            plt.close()

            graph_nshare(all_data_series)
            pdf.savefig()
            plt.close()

    else:
        graph_gc(all_data_series)
        plt.show()
        graph_len_count(all_data_series)
        plt.show()
        graph_scores(all_data_series)
        plt.show()
        graph_basequality(all_data_series)
        plt.show()
        graph_basenvert(all_data_series)
        plt.show()
        graph_basenvert_abr(all_data_series)
        plt.show()
        graph_nshare(all_data_series)
        plt.show()

    tabelle_speichern(all_data_series, arguments.save_table)
    print("Table saved.")

if __name__ == "__main__":
    start = time.time()
    main()
    end = time.time()
    print(end-start)
