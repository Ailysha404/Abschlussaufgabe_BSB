#unterscheidung phred 64 or 32
#dictionary mit alphabet zu qualität
#aber nich gesamte sequenz überprüfen, ab gewissem unsicherheitspunkt fragen
#qualitätskontrolle, cut-off vom user angegeben, sonst rausgeworfen
#bei trimming bestimmen härte des cut-offs
#evtl erst ab nachfolgenden schlechten qualitäten, oder mittelwert fenster
#auswertung für: verteilung der readlängen mit mw, median, stabw
#gc-gehalt
#qualitätsdaten mit mw, pro base, pro sequenz, gesamtfenster pro base
#graphische darstellung, evtl plotly?
#evtl prüfen auf Ns, evtl k-mere?
#optionsmanagement über kommandozeilen-parameter - interaktiver und nicht interaktiver modus - bei standard KEINE nutzerabfragen
#evtl datei mit output speichern zur nachverfolgung
#präsentation 21.02, ~20 min als team, abgabe am tag selbst

import sys
import argparse
import matplotlib.pyplot as plt

def Phred_Bestimmung(datei):
    alphabet = {} #Alphabet ascii-zeichen, score
    score_raw = "" #read line, string ascii
    pass

def Parser():
    input_parser = argparse.ArgumentParser(description="Qualitätsauswertung für FASTA-Dateien")
    input_parser.add_argument("Dateipfad", help="Dateipfad zur FASTA-Datei")
    input_parser.add_argument('-p','--phred', dest='opt', help="Optionale manuelle Angabe des Phred-Formats, 33 oder 64")
    input_parser.add_argument('-i','--interaktiv', dest='interaktiv', action="store_true", help="Aktivierung des kommentierten interaktiven Modus")
    arguments = input_parser.parse_args() #Erstellen dictionary, Zugriff auf Argument über Namen
    return arguments
    pass

def Qualitaet(datei,phred=None):
    values = {}
    ascii = ("!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
    alphabet_33 = {char:value-33 for value,char in enumerate(ascii,33) if value < 74} #make dictionary w ascii-character + phred score
    alphabet_64 = {char:value-64 for value,char in enumerate(ascii,33) if value > 63 and value < 105} #make dictionary w ascii-character + phred score
    scores_ascii = dict(zip(datei[0::4],datei[3::4])) #keys ids, values score
    if phred == None: #if not manually given phred, check first sequence for clarification
        for score_char in scores_ascii.values():
            if score_char in alphabet_33.keys():
                phred = "33"

            elif score_char in alphabet_64.keys():
                phred = "64"

            else:
                quit()

    if phred == "33": #check for phred format
        for score_string in scores_ascii: #id sequence
            for letter in scores_ascii[score_string]: #individual letter in score line
                print(f"{letter}, {alphabet_33[str(letter)]}")

    elif phred == "64":
        for score_string in scores_ascii:
            score_list = []
            for letter in scores_ascii[score_string].rstrip(): #return score for each ascii
                score_list.append(alphabet_64[str(letter)]) #
                values[score_string+letter] = score_list

    return values

def __main__():

    arguments = Parser()
    file = open(arguments.Dateipfad)
    inhalt = file.readlines()[0:8] #lines als Liste gespeichert - 0:id 1:seq 2:id 3:score und multiple...
    file.close()
    print(Qualitaet(inhalt,arguments.opt))
    #scores_seq = Qualitaet(inhalt) #Dictionary mit key ids, values score
    #print(f"{scores_seq.keys()} : {scores_seq.values()}")
    pass

if __name__ == "__main__":
    __main__()
