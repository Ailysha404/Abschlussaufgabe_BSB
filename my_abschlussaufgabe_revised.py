import argparse
import matplotlib.pyplot as plt
import time
import numpy as np
from tabulate import tabulate
import pandas as pd
import seaborn as sns
import csv
from itertools import zip_longest
import matplotlib.ticker as ticker


class Read:
    """Store read with all necessary information as attributes."""
    def __init__(self, content, id, phred, alphabet, quality):
        self.length = len(content[1])
        self.gc = (content[1].count("G") + content[1].count("C"))/self.length #in decimal - between 0 and 1
        self.name = content[0].split()[0] #Id from first line as name
        self.quality = quality
        self.sequenz = content[1][0:len(self.quality)] #sequence-string


def Parser():
    """Parses command line input when starting script."""
    input_parser = argparse.ArgumentParser(description="Qualitätsauswertung für FASTA-Dateien")
    input_parser.add_argument("Dateipfad", help="Dateipfad zur FASTA-Datei")
    input_parser.add_argument('-p','--phred', dest='phred', default=None, help="Optionale manuelle Angabe des Phred-Formats, 33 oder 64")
    input_parser.add_argument('-t','--trim', dest='trim_val', help="Optionale manuelle Angabe des Trimming-Scores",default=25)
    input_parser.add_argument('-c','--cutoff', dest='cutoff', help="Optionale manuelle Angabe des Durchschnittscores, bei dem ein Read komplett verworfen wird.",default=20)
    input_parser.add_argument('-m','--minl', dest='minlength', help="Optionale manuelle Angabe der minimalen Länge in Basen, die ein getrimmter Read bei der Auswertung haben muss.",default=20)
    input_parser.add_argument('-tab','--table', dest='save_table', help="Optionale Angabe eines Dateipfades zum Speichern der Datentabelle", default=None)
    input_parser.add_argument('-plt','--plot', dest='save_plot', help="Optionale Angabe eines Dateipfades zum Speichern der Datenplots", default=None)
    input_parser.add_argument('-i','--interaktiv', dest='interaktiv', action="store_true", help="Aktivierung des kommentierten interaktiven Modus")
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


def Phred_Bestimmung(datei, ascii, phred=None):
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


def Qualitaet(content, alphabet, phred):
    """Create Phred-Score dictionary and translate ASCII-character into score.
    Append score to value list.
    """
    values = []
    values += [alphabet[str(letter)] for letter in content[3]]

    return values


def trimming(scores,trim_val):
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
            summe = 0 #Summe wieder zurücksetzen!

            if kmer_mean < trim_val:
                return scores[0:item_index]

    else:
        return scores


def graph_gc(dict_reads, fig):
    """Create subplot for gc-percentages per read."""
    fig.add_subplot(4, 2, (1, 2))
    ax = sns.countplot([item.gc*100 for item in dict_reads])
    ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set(xlabel="GC content in %",ylabel="Number Reads",title="GC content")

    return ax


def graph_len_count(dict_reads, fig):
    fig.add_subplot(4, 2, (3, 4))
    ax = sns.countplot([item.length for item in dict_reads], palette="pastel")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set(xlabel="Readlength in BP",ylabel="Number Reads",title="Readlengths")

    return ax


def graph_scores(dict_reads, fig):
    """Create subplot for mean scores per read."""
    fig.add_subplot(4, 2, (5, 6))
    x_val = [round(sum(score*10 for score in item.quality)/len(item.quality))/10 for item in dict_reads]
    ax = sns.countplot(x_val)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set(xlabel="Quality Score",ylabel="Number Reads",title="Sequence Quality")

    return ax


def graph_basequality(dict_reads, fig):
    """Shows the Quality of all Read positions"""
    fig.add_subplot(4, 2, (7, 8))
    qualscores = [item.quality for item in dict_reads]
    array_qual = np.array(list(zip_longest(*qualscores,fillvalue=np.nan)), dtype=float)
    dataframe = pd.DataFrame(array_qual.tolist())

    #FERTIG!!SIEEEG!!
    sterr = np.std(dataframe, axis=1)
    ax = plt.errorbar(x=dataframe.index, y=[np.mean(dataframe.loc[i]) for i in dataframe.index], yerr=sterr, data=dataframe, ecolor="cyan", capthick="0.7")

    return ax


def graph_basenvert(dict_reads, fig):
    """Shows the Quality of all Read positions"""
    fig.add_subplot(4, 2, (1, 2))
    sequences = [item.sequenz for item in dict_reads]
    dataframe = list(zip_longest(*sequences,fillvalue="X"))

    ax = sns.lineplot(x=[index for index in range(0,max([len(list(dataframe))]))], y=list([i.count("A")/len(''.join(i).replace("X","").strip(", ")) for i in list(dataframe)]))
    ax = sns.lineplot(x=[index for index in range(0,max([len(list(dataframe))]))], y=list([i.count("C")/len(''.join(i).replace("X","").strip(", ")) for i in list(dataframe)]))
    ax = sns.lineplot(x=[index for index in range(0,max([len(list(dataframe))]))], y=list([i.count("G")/len(''.join(i).replace("X","").strip(", ")) for i in list(dataframe)]))
    ax = sns.lineplot(x=[index for index in range(0,max([len(list(dataframe))]))], y=list([i.count("T")/len(''.join(i).replace("X","").strip(", ")) for i in list(dataframe)]))
    ax = sns.lineplot(x=[index for index in range(0,max([len(list(dataframe))]))], y=list([i.count("N")/len(''.join(i).replace("X","").strip(", ")) for i in list(dataframe)]))

    return ax


def tabelle_speichern(dict_reads, dateipfad="Datentabelle.csv"):
    # '''Tabelle speichern in csv'''
    dataframe = pd.DataFrame({"Name":item.name, "Länge":item.length, "GC":item.gc, "Qualität":str(item.quality).strip("[]"), "Sequenz":item.sequenz} for item in dict_reads)
    dataframe.index = np.arange(1,len(dataframe)+1)
    f = open("table.csv","w")
    dataframe.to_csv(path_or_buf=dateipfad,header=dataframe.columns,index_label="Index")


def main():
    all_ids = [] #! saved as object code, not as string!
    ascii = ("!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
    arguments = Parser()
    cutoff = arguments.cutoff
    minlength = arguments.minlength
    phred = arguments.phred
    trim_val = arguments.trim_val

    with open(arguments.Dateipfad) as inhalt:
        block = []
        test_cnt = 0	#!! KÜRZEN FÜR TESTLÄUFE

        for lines in inhalt:

            if block == []:
                phred = Phred_Bestimmung(lines, ascii, phred)

            block.append(lines.rstrip())
            test_cnt += 1	#!! KÜRZEN FÜR TESTLÄUFE

            if test_cnt > 1000:	#!! KÜRZEN FÜR TESTLÄUFE
                break

    alphabet = alphabet_ascii(phred, ascii)
    line_pack = []

    for lines in block:
        line_pack.append(lines)

        if len(line_pack) == 4:
            quality = trimming(Qualitaet(line_pack, alphabet, phred),trim_val)
            if np.mean(quality) > cutoff and len(quality) >= minlength:
                all_ids.append(Read(line_pack, id, phred, alphabet, quality))
            line_pack = []

    #Following code is responsible for creation and saving of graphs
    fig = plt.figure()
    sns.set_style("whitegrid")

    # graph_gc(all_ids, fig)
    # graph_scores(all_ids, fig)
    # graph_len_count(all_ids, fig)
    # graph_basequality(all_ids, fig)
    # graph_basenvert(all_ids, fig)

    #gc_boxplot(all_ids, fig)
    #
    # if arguments.save_plot:
    # 	plt.savefig(str(arguments.save_plot))
    # else:
    #     plt.show()

    tabelle_speichern(all_ids, arguments.save_table)

if __name__ == "__main__":
    #start = time.time()
    main()
    #end = time.time()

# https://pypi.org/project/tabulate/
# print(end-start)
# test_table = {"c":[500, 500, 30], "b":[30, 40, 2], "a":[40, 40, 5]}
# print(tabulate((test_table), headers="keys"))
# file = open("/home/ailysha/abschluss/table.txt","w+")
# file.write(tabulate((test_table), headers="keys"))

#seaborn
# sns.set(style="darkgrid")
# set = dataset
# sns.relplot(x="timepoint", y="signal", kind="line", data=set)
# plt.show()
