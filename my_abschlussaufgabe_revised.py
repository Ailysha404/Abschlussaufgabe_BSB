import argparse
import matplotlib.pyplot as plt


class Read:
    """Store read with all necessary information as attributes."""
    def __init__(self, content, id, phred, alphabet):
        self.length = len(content[1])
        self.gc = (content[1].count("G") + content[1].count("C"))/self.length #in decimal - between 0 and 1
        self.name = content[0].split()[0] #Id from first line as name
        self.quality = trimming(self.qualitaet(content, alphabet, phred))
        self.sequenz = content[1][0:len(self.quality)] #sequence-string

    def qualitaet(self, content, alphabet, phred):
        values = []

        if phred == "33": #check for phred format
            alphabet_33 = {char:value for value, char in enumerate(alphabet)} #make dictionary w ascii-character + phred score

            for letter in alphabet_33[score_string].rstrip(): #individual letter in score line
                    values.append(alphabet_33[str(letter)]) #add score to list

        elif phred == "64":
            alphabet_64 = {char:value for value, char in enumerate(alphabet)} #make dictionary w ascii-character + phred score

            for letter in content[3]:
                    values.append(alphabet_64[str(letter)])

        else:
            quit()

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
    alphabet = None

    if phred == "33": #check for phred format
        alphabet = {char:value-33 for value,char in enumerate(ascii,33) if value < 74} #make dictionary w ascii-character + phred score

    elif phred == "64":
        alphabet = {char:value-64 for value,char in enumerate(ascii,33) if value > 63 and value < 105} #make dictionary w ascii-character + phred score

    return alphabet


def Parser():
    input_parser = argparse.ArgumentParser(description="Qualitätsauswertung für FASTA-Dateien")
    input_parser.add_argument("Dateipfad", help="Dateipfad zur FASTA-Datei")
    input_parser.add_argument('-p','--phred', dest='phred', default=None, help="Optionale manuelle Angabe des Phred-Formats, 33 oder 64")
    input_parser.add_argument('-t','--trim', dest='trim_val', help="Optionale manuelle Angabe des Trimming-Scores",default=25)
    input_parser.add_argument('-c','--cutoff', dest='cutoff', help="Optionale manuelle Angabe des Durchschnittscores, bei dem ein Read komplett verworfen wird.",default=28)
    input_parser.add_argument('-m','--minl', dest='minlength', help="Optionale manuelle Angabe der minimalen Länge in Basen, die ein getrimmter Read bei der Auswertung haben muss.",default=20)
    input_parser.add_argument('-s','--save', dest='save', help="Optionale Angabe eines Dateipfades zum Speichern der Daten")
    input_parser.add_argument('-i','--interaktiv', dest='interaktiv', action="store_true", help="Aktivierung des kommentierten interaktiven Modus")
    arguments = input_parser.parse_args() #Erstellen dictionary, Zugriff auf Argument über Namen

    return arguments


def trimming(scores):
    counter = 0
    summe = 0
    item_index = 0
    kmer_neu = []

    for kmer in scores:
        counter += 1
        summe += kmer
        item_index += 1

        if counter == 3:
            kmer_mean = summe / 3
            counter = 0
            summe = 0 #Summe wieder zurücksetzen!

            if kmer_mean < 25:
                kmer_neu = scores[0:item_index]
                return kmer_neu

    else:
        return scores

def graph_gc(dict_reads, fig):
    fig.add_subplot(3, 2, (1, 2))
    x_values = []

    for value in dict_reads:
    	x_values.append(value.gc * 100)

    plt.hist(x_values, histtype="step", color="blue")
    plt.xlabel("GC Content (%)")
    plt.ylabel("Count")
    plt.title("GC Distribution over all Sequences")

    return plt

def graph_scores(dict_reads, fig):
    fig.add_subplot(3, 2, (5, 6))
    x_values = []

    for item in dict_reads:
        overall_score = 0

        for score in item.quality:
            overall_score += score

        x_values.append(overall_score / len(item.quality))

    plt.hist(x_values,histtype="step",color="orange",label="Average Quality per Read",lw=2)
    plt.xlabel("Mean Sequence Quality")
    plt.ylabel("Reads")
    plt.title("Quality Score Distribution over all Sequences")

    return plt

def main():
    all_ids = [] #! saved as object code, not as string! saved in item.name
    ascii = ("!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
    arguments = Parser()
    phred = arguments.phred
    block = []

    with open(arguments.Dateipfad) as inhalt: #open user-provided file path
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
            id = line_pack[0]
            all_ids.append(Read(line_pack, id, phred, alphabet))
            line_pack = []

    fig = plt.figure()	#create canvas for graphs
    graph_gc(all_ids, fig)	#create first subplot, for gc-percentages
    graph_scores(all_ids, fig)	#fill in second subplot, for scores
    plt.show()
    # if arguments.save:
    # 	plt.savefig(str(arguments.save))

if __name__ == "__main__":
    main()
