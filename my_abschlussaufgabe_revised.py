import argparse
import matplotlib.pyplot as plt
import time
import numpy as np
from tabulate import tabulate
import pandas as pd
import seaborn as sns


class Read:
    """Store read with all necessary information as attributes."""
    def __init__(self, content, id, phred, alphabet, quality):
        self.length = len(content[1])
        self.gc = (content[1].count("G") + content[1].count("C"))/self.length #in decimal - between 0 and 1
        self.name = content[0].split()[0] #Id from first line as name
        self.quality = quality
        self.sequenz = content[1][0:len(self.quality)] #sequence-string



def qualitaet(content, alphabet, phred):
    """Create Phred-Score dictionary and translate ASCII-character into score.
    Append score to value list.
    """
    values = []
    values += [alphabet[str(letter)] for letter in content[3]]

    return values


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


def alphabet_ascii(phred, ascii):
    """Create dictionary with scores for ASCII by given phred format."""
    alphabet = None

    if phred == "33": #check for phred format
        alphabet = {char:value-33 for value,char in enumerate(ascii,33) if value < 74} #make dictionary w ascii-character + phred score

    elif phred == "64":
        alphabet = {char:value-64 for value,char in enumerate(ascii,33) if value > 63 and value < 105} #make dictionary w ascii-character + phred score

    return alphabet


def Parser():
    """Parses command line input when starting script."""
    input_parser = argparse.ArgumentParser(description="Qualitätsauswertung für FASTA-Dateien")
    input_parser.add_argument("Dateipfad", help="Dateipfad zur FASTA-Datei")
    input_parser.add_argument('-p','--phred', dest='phred', default=None, help="Optionale manuelle Angabe des Phred-Formats, 33 oder 64")
    input_parser.add_argument('-t','--trim', dest='trim_val', help="Optionale manuelle Angabe des Trimming-Scores",default=25)
    input_parser.add_argument('-c','--cutoff', dest='cutoff', help="Optionale manuelle Angabe des Durchschnittscores, bei dem ein Read komplett verworfen wird.",default=20)
    input_parser.add_argument('-m','--minl', dest='minlength', help="Optionale manuelle Angabe der minimalen Länge in Basen, die ein getrimmter Read bei der Auswertung haben muss.",default=20)
    input_parser.add_argument('-s','--save', dest='save', help="Optionale Angabe eines Dateipfades zum Speichern der Daten")
    input_parser.add_argument('-i','--interaktiv', dest='interaktiv', action="store_true", help="Aktivierung des kommentierten interaktiven Modus")
    arguments = input_parser.parse_args() #Erstellen dictionary, Zugriff auf Argument über Namen

    return arguments


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
    fig.add_subplot(3, 2, (1, 2))
    x_values = []

    x_values = [item.gc*100 for item in dict_reads]
    ax = sns.distplot(x_values,kde=False,color="orange")
    ax.set(xlabel="GC content in %",ylabel="Number Reads",title="GC content")

    return ax

def graph_scores(dict_reads, fig):
    """Create subplot for mean scores per read."""
    fig.add_subplot(3, 2, (3, 4))
    x_values = []

    x_values = [sum(score for score in item.quality)/len(item.quality) for item in dict_reads]
    ax = sns.distplot(x_values,kde=False,color="orange")
    ax.set(xlabel="Quality Score",ylabel="Number Reads",title="Sequence Quality")

    return ax

def graph_basequality(dict_reads, fig, dataframe):
    """Shows the Quality of all Read positions"""
    fig.add_subplot(3, 2, (5, 6))
    #Ich hasse Seaborn und pandas und matplotlib


    #X=[score_list]#x and y have to be the same lenght
    #Y=[range(0,len(score_list))]
    #y_values = [item.info["Qualität"] for item in dict_reads]
    #x_values = [index for index, item in enumerate(y_values)]

    #ax=sns.catplot("Name", "Länge", col="GC", data=dataframe, kind="point")

    #ax = plt.errorbar(dataframe.index, dataframe["GC"], yerr=20 , data=dataframe)
    # print(len(range(0,max(dataframe["Länge"])+1)))
    # print(sum(dataframe["Qualität"][1])/len(dataframe["Qualität"][1]))


    #print(len([sum(dataframe["Qualität"][i]) / len(dataframe["Qualität"][i]) for i in range(0,max(dataframe["Länge"])+1)]))


    # ax = sns.lineplot(x=range(0,max(dataframe["Länge"])+1), y=[sum(dataframe["Qualität"][i]) / len(dataframe["Qualität"][i]) for i in range(0,max(dataframe["Länge"]+1))], data=dataframe)


    # plt.plot(Y,X,alpha=0.7,color="black")
    # plt.xlabel('Posittion')
    # plt.ylabel('Quality')
    # plt.title('Quality Distribution over all Reads')

    return ax

def main():
    all_ids = [] #! saved as object code, not as string!
    ascii = ("!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
    arguments = Parser()
    cutoff = arguments.cutoff
    minlength = arguments.minlength
    phred = arguments.phred
    trim_val = arguments.trim_val
    block = []

    with open(arguments.Dateipfad) as inhalt:
        test_cnt = 0	#!! KÜRZEN FÜR TESTLÄUFE

        for lines in inhalt:

            if block == []:
                phred = Phred_Bestimmung(lines, ascii, phred)

            block.append(lines.rstrip())
            test_cnt += 1	#!! KÜRZEN FÜR TESTLÄUFE

            if test_cnt > 500:	#!! KÜRZEN FÜR TESTLÄUFE
                break

    alphabet = alphabet_ascii(phred, ascii)
    line_pack = []

    for lines in block:
        line_pack.append(lines)

        if len(line_pack) == 4:
            quality = trimming(qualitaet(line_pack, alphabet, phred),trim_val)
            if np.mean(quality) > cutoff and len(quality) >= minlength:
                all_ids.append(Read(line_pack, id, phred, alphabet, quality))
            line_pack = []

    test = ({"Länge":item.length, "GC":item.gc, "Name":item.name, "Qualität":item.quality, "Sequenz":item.sequenz} for item in all_ids)
    info = pd.DataFrame(test)

    #Following code is responsible for creation and saving of graphs
    fig = plt.figure()
    sns.set_style("whitegrid")
    graph_gc(all_ids, fig)
    graph_scores(all_ids, fig)
    graph_basequality(all_ids, fig, info)
    plt.show()
    if arguments.save:
    	plt.savefig(str(arguments.save))

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
