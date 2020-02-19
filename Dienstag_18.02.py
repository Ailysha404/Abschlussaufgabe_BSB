import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
import time
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
        default=None
        )
    input_parser.add_argument(
        "-plt", "--plot",
        dest="save_plot",
        help="Optionale Angabe eines Dateipfades zum Speichern der Datenplots",
        default=None
        )
    input_parser.add_argument(
        "-i", "--interaktiv",
        dest="interaktiv",
        action="store_true",
        help="Aktivierung des kommentierten interaktiven Modus"
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


def graph_gc(dict_reads, fig):
    """Create subplot for gc-percentages per read."""
    # gs = GridSpec(1, 4, figure = fig)
    # pdata = pd.DataFrame(data=pd.Index([int(round(item.gc,2)*100) for item in dict_reads]).value_counts())
    #
    # ax = fig.add_subplot(gs[0,:], label="Readlength in BP")
    # sns.lineplot(data=pdata, legend=None, palette=["teal"])

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
        )

    return fig


def graph_len_count(dict_reads, fig):
    """Create graph showing read length distribution over the reads."""
    maxLen = max([item.length for item in dict_reads])
    # gs = GridSpec(2, 2, figure = fig, wspace=.5, hspace=.5)
    # pdata = pd.DataFrame(data=pd.Index([item.length for item in dict_reads]).value_counts()).sort_index().sub(1,axis=1)

    ax1, ax2 = fig.subplots(2, 1, sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.2)

    ax1.set(
        title="Readlengths before Trimming",
        ylabel="Number of Reads",
        )
    ax1.set_xlim(0, maxLen)
    sns.distplot(
        [item.length for item in dict_reads],
        bins=40,
        kde=False,
        ax=ax1,
        color="teal",
        )

    ax2.set(
        xlabel="Readlength in BP",
        ylabel="Number of Reads",
        )

    sns.distplot(
        [item.trimmed_length for item in dict_reads],
        bins=40,
        kde=False,
        ax=ax2,
        color="teal",
        )

    #sns.lineplot(data=pdata, legend=None, palette=["teal"])
    return fig


def graph_scores(dict_reads, fig):
    """Create subplot for mean scores per read."""
    # gs = GridSpec(1, 4, figure = fig)
    # ax = fig.add_subplot(gs[0,:])
    # x_val = [round(sum(score*10 for score in item.quality)/len(item.quality))/10 for item in dict_reads]
    #
    # ax = sns.countplot(x_val, palette="viridis_r",)
    # ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=10))
    # plt.xlabel("Quality Score")
    # plt.ylabel("Number of Reads")
    # plt.title("Sequence Quality")

    ax1, ax2 = fig.subplots(2, 1, sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.2)

    ax1.set(
        title="Average Sequence Quality",
        ylabel="Number of Reads",
        )
    sns.distplot(
    [item.mean_qual_untrimmed for item in dict_reads],
        bins=40,
        kde=False,
        ax=ax1,
        color="teal",
        )

    ax2.set(
        xlabel="Quality Score",
        ylabel="Number of Reads",
        )
    sns.distplot(
    [item.mean_qual_trimmed for item in dict_reads],
        bins=40,
        kde=False,
        ax=ax2,
        color="teal",
        )

    return fig


def graph_basequality(dict_reads, fig):
    """Show the quality of all Read positions"""
    # gs = GridSpec(1, 4, figure = fig)
    # fig.add_subplot(gs[0,:])
    qualscores = [item.quality for item in dict_reads]
    array_qual = np.array(list(zip_longest(*qualscores,fillvalue=np.nan)), dtype=float)
    dataframe = pd.DataFrame(array_qual.tolist())
    sterr = np.std(dataframe, axis=1)
    ax = fig.subplots(1, 1, sharex=True, sharey=True)
    ax.set(
        title="Average Quality at Read Positions",
        xlabel="Position in Read (BP)",
        ylabel="Quality",
        )
    plt.errorbar(
        x=dataframe.index,
        y=[np.mean(dataframe.loc[i]) for i in dataframe.index],
        yerr=sterr,
        data=dataframe,
        color="#008080",
        lw=0.8,
        ecolor="lightblue",
        capthick="0.5",
        )

    return fig


def graph_basenvert(dict_reads, fig):
    """Show the distribution of all bases"""
    # # gs = GridSpec(1, 4, figure = fig)
    # sns.set_style({"palette":"Dark2"})
    # # fig.add_subplot(gs[0,:])
    # sequences = [item.sequenz for item in dict_reads]
    # datalist = list(zip_longest(*sequences,fillvalue=np.nan))
    # colors=["red", "green", "blue", "yellow", "purple"]
    # bases = ["A", "C", "G", "T", "N"]
    # y_values = []
    # #x_values =
    #
    # counter = 0
    #
    # ax = fig.subplots(1, 1, sharex=True, sharey=True)
    # ax.set(
    #     title="Distribution of Bases at Read Positions",
    #     xlabel="Position in Read (BP)",
    #     ylabel="Distribution of Bases (%)",
    #     )
    #
    # for letter in bases:
    #     for group in [i.count(letter)/len(list(x for x in i if not pd.isna(x)))*100 for i in datalist]:
    #         y_values.append(np.nanmean(group))
    #     letter = sns.lineplot(
    #         y=y_values,
    #         x=[index for index in range(0,len(y_values))],
    #         lw=0.5,
    #         alpha=0.5,
    #         color=colors[counter]
    #         )
    #     counter += 1
    #
    #
    # fig.legend(
    # bases,
    # loc="upper center",
    # ncol=4,
    # fancybox=True,
    # framealpha=1,
    # edgecolor="#000000",
    # facecolor="#FFFFFF",
    # )

    sequences = [list(item.sequenz) for item in dict_reads]
    datalist = list(zip_longest(*sequences, fillvalue=np.nan))
    bases = ["A", "C", "G", "T", "N"]
    colors=["red", "green", "blue", "yellow", "purple"]
    counter = 0
    means = []
    calc = []
    #ax = fig.subplots(1, 1)
    ax = fig.subplots(1, 1, sharex=True, sharey=True)
    ax.set(
        title="Distribution of Bases at Read Positions",
        xlabel="Position in Read (BP)",
        ylabel="Distribution of Bases (%)",
        )
    for letter in bases:
        y_values = []
        for group in [i.count(letter)/len(list(x for x in i if not pd.isna(x)))*100 for i in datalist]:
            calc.append(group)
            #print(list(chain(calc)))
            if len(calc) > 4:
                y_values.append(np.nanmean(list(chain(calc))))
                calc=[]

        letter = sns.lineplot(
            x=[index*4 for index in range(0,max([len(y_values)]))],
            y=y_values,
            lw=0.5,
            color=colors[counter],
            legend=False,
            )
        counter += 1
    return fig


def graph_basenvert_abr(dict_reads, fig, pdf):
    """Show the distribution of all bases"""
    sequences = [list(item.sequenz) for item in dict_reads]
    datalist = list(zip_longest(*sequences, fillvalue=np.nan))
    bases = ["A", "C", "G", "T", "N"]
    colors=["red", "green", "blue", "yellow", "purple"]
    counter = 0
    means = []
    calc = []

    for letter in bases:
        ax = fig.subplots(1, 1)
        y_values = []
        for group in [i.count(letter)/len(list(x for x in i if not pd.isna(x)))*100 for i in datalist]:
            calc.append(group)
            #print(list(chain(calc)))
            if len(calc) > 4:
                y_values.append(np.nanmean(list(chain(calc))))
                calc=[]

        plt.title(f"Average Percentage of {letter} at Read Positions")
        plt.xlabel("Position in Read (BP)")
        plt.ylabel(f"Share of {letter} over all Reads (%)")
        plt.ylim(0, 100)
        ax.set_xticklabels
        letter = sns.lineplot(
            x=[index*4 for index in range(0,max([len(y_values)]))],
            y=y_values,
            lw=0.5,
            color=colors[counter],
            legend=False,
            )
        counter += 1
        pdf.savefig()
        plt.close()
    # sequences = [item.sequenz for item in dict_reads]
    # datalist = list(zip_longest(*sequences,fillvalue="X"))
    # x_values = [index for index in range(0,max([len(datalist)]))]
    # y_values = []
    # bases = ["A", "C", "G", "T"]
    # count = 0
    # gs = GridSpec(4, 1, figure = fig)
    # fig.subplots(4, 1, sharex=True, sharey=True)
    # fig.subplots_adjust(hspace=0.5, wspace=0.2)
    # for letter in bases:
    #
    #     ax = fig.add_subplot(gs[count,:])
    #     ax.set(
    #         title=f"Percentage of {letter} at Read Positions",
    #         xlabel="Position in Read (BP)",
    #         ylabel=f"Share of {letter} over all Reads (%)",
    #         )
    #     divisor=5
    #     gesamtlänge= [i.count(letter) for i in datalist]
    #     x_val = [round(sum(score*10 for score in item.quality)/len(item.quality))/10 for item in dict_reads]
    #     y_values = list([i.count(letter)/len("".join(i).replace("X","").strip(", "))*100 for i in datalist])
    #     count += 1
    #     #
    #     # ax = sns.countplot(x_val, palette="viridis_r",)
    #     # ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=10))
    #     # for item in range(1,10):
    #     #     if item % len(gesamtlänge) == 0:
    #     #         divisor = item
    #     #         y_values = dataframe[current:current + divisor]/divisor for current in range(0,gesamtlänge,divisor)
    #     print(y_values)
    #     print(x_values)
    #     letter = sns.lineplot(
    #         x=x_values,
    #         y=y_values,
    #         alpha=0.7,
    #         lw=0.5,
    #         ax=ax,
    #         )
    #     pdf.savefig()
    #
    # fig.legend(
    # bases,
    # loc="upper center",
    # ncol=4,
    # fancybox=True,
    # framealpha=1,
    # edgecolor="#000000",
    # facecolor="#FFFFFF",
    # )

    return fig


def tabelle_speichern(dict_reads, phred, read_count, dateipfad="Datentabelle.csv"):
    """Save data stats/attributes as csv table for later review"""
    dataframe_overview = pd.DataFrame({
        "Anzahl Reads":[
            read_count,
            "",
            len(dict_reads),
            ""
            ],
        "Phred-Score Sytem":phred,
        "Durchschn. Länge":[
            np.mean([item.length for item in dict_reads]),
            np.std([item.length for item in dict_reads]),
            np.mean([item.trimmed_length for item in dict_reads]),
            np.std([item.trimmed_length for item in dict_reads])
            ],
        "Durchschn. Qualität":[
            np.mean([item.mean_qual_untrimmed for item in dict_reads]),
            np.std([item.mean_qual_untrimmed for item in dict_reads]),
            np.mean([item.mean_qual_trimmed for item in dict_reads]),
            np.std([item.mean_qual_trimmed for item in dict_reads])
            ],
        "Durchschn. GC-Gehalt %":[
            "",
            "",
            np.mean([item.gc*100 for item in dict_reads]),
            np.std([item.gc*100 for item in dict_reads])
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
    all_ids = [] #! saved as object code, not as string!
    ascii = ("!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
    arguments = parser()
    cutoff = arguments.cutoff
    minlength = arguments.minlength
    phred = arguments.phred
    trim_val = arguments.trim_val

    with open(arguments.Dateipfad) as inhalt:
        block = []
        test_cnt = 0	#!! KÜRZEN FÜR TESTLÄUFE

        for lines in inhalt:

            if block == []:
                phred = phred_bestimmung(
                    lines,
                    ascii,
                    phred
                    )

            block.append(lines.rstrip())
            test_cnt += 1	#!! KÜRZEN FÜR TESTLÄUFE

            if test_cnt > 1000:	#!! KÜRZEN FÜR TESTLÄUFE
                break

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

    #Following code is responsible for creation and saving of graphs
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

            fig = plt.figure()
            graph_gc(all_ids, fig)
            pdf.savefig()
            plt.close()

            fig = plt.figure()
            fig = graph_len_count(all_ids, fig)
            pdf.savefig()
            plt.close()

            fig = plt.figure()
            graph_scores(all_ids, fig)
            pdf.savefig()
            plt.close()

            fig = plt.figure()
            graph_basequality(all_ids, fig)
            pdf.savefig()
            plt.close()

            fig = plt.figure()
            graph_basenvert(all_ids, fig)
            pdf.savefig()
            plt.close()

            fig = plt.figure()
            graph_basenvert_abr(all_ids, fig, pdf)



    else:
        fig = plt.figure()
        graph_gc(all_ids, fig)
        plt.show()
        graph_len_count(all_ids, fig)
        plt.show()
        graph_scores(all_ids, fig)
        plt.show()
        graph_basequality(all_ids, fig)
        plt.show()
        graph_basenvert(all_ids, fig)
        plt.show()
        graph_basenvert_abr(all_ids, fig)
        plt.show()

    tabelle_speichern(all_ids, phred, read_count, arguments.save_table)

if __name__ == "__main__":
    start = time.time()
    main()
    end = time.time()
    print(end-start)



# https://pypi.org/project/tabulate/
# test_table = {"c":[500, 500, 30], "b":[30, 40, 2], "a":[40, 40, 5]}
# print(tabulate((test_table), headers="keys"))
# file = open("/home/ailysha/abschluss/table.txt","w+")
# file.write(tabulate((test_table), headers="keys"))

#seaborn
# sns.set(style="darkgrid")
# set = dataset
# sns.relplot(x="timepoint", y="signal", kind="line", data=set)
# plt.show()
