import sys
import argparse
import matplotlib.pyplot as plt

class Read: #class for all reads
	def __init__(self,content,id,phred,ascii,idlist):
		self.length = len(content[1])
		self.gc = (content[1].count("G")+content[1].count("C"))/self.length #in decimal - between 0 and 1
		self.name = content[0].split()[0] #how to get specific substring from line 1? #saves own id as name

		self.quality = trimming(self.qualitaet(content,phred,ascii)[self]) #dict of name with quality-scores
		self.sequenz = content[1][0:len(self.quality)] #sequence-string

	def qualitaet(self,content,ascii,phred):
		values = {}

		if phred == "33": #check for phred format
			alphabet_33 = {char:value-33 for value,char in enumerate(ascii,33) if value < 74} #make dictionary w ascii-character + phred score

			for score_string in alphabet_33: #id read
				score_list = [] #list with all scores from the read

				for letter in alphabet_33[score_string].rstrip(): #individual letter in score line
					score_list.append(alphabet_33[str(letter)]) #add score to list
					values[self] = score_list  #dictionary with id string as key, list with scores as value

		elif phred == "64":
			alphabet_64 = {char:value-64 for value,char in enumerate(ascii,33) if value > 63 and value < 105} #make dictionary w ascii-character + phred score

			for score_string in alphabet_64:
				score_list = []

				for letter in content[3]:
					score_list.append(alphabet_64[str(letter)])
					values[self] = score_list

		else:
			quit()

		return values

def Phred_Bestimmung(datei,ascii,phred=None):

	if phred == None: #if not manually given phred, check first sequence for clarification
		for score_char in datei:

			if score_char in [char for value,char in enumerate(ascii,33) if value < 64]:
				phred = "33"
				break

			elif score_char in [char for value,char in enumerate(ascii,33) if value > 63 and value < 105]:
				phred = "64"
				break

			else:
				quit()

	return phred

def Parser():
	input_parser = argparse.ArgumentParser(description="Qualitätsauswertung für FASTA-Dateien")
	input_parser.add_argument("Dateipfad", help="Dateipfad zur FASTA-Datei")
	input_parser.add_argument('-p','--phred', dest='phred', help="Optionale manuelle Angabe des Phred-Formats, 33 oder 64")
	input_parser.add_argument('-t','--trim', dest='trim_val', help="Optionale manuelle Angabe des Trimming-Scores",default=25)
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
			kmer_mean=summe/3
			counter=0
			summe=0 #Summe wieder zurücksetzen!

			if kmer_mean < 25:
				kmer_neu = scores[0:item_index]
				return kmer_neu
	else:
		return scores

def graph_gc(dict_reads,fig):
	fig.add_subplot(4,2,(1,4))
	X=[]
	for value in dict_reads:
		X.append(value.gc*100)

	plt.hist(X,histtype="step",color="blue")
	plt.xlabel('GC Content (%)')
	plt.ylabel('Count')
	plt.title('GC Distribution over all Sequences')
	return plt

def graph_scores(dict_reads,fig):
	fig.add_subplot(4,2,(5,8))
	X=[]
	for item in dict_reads:
		overall_score = 0
		for score in item.quality:
			overall_score += score
		X.append(overall_score/len(item.quality))

	plt.hist(X,histtype="step",color="orange",label="Average Quality per Read",lw=2)
	plt.xlabel('Mean Sequence Quality')
	plt.ylabel('Reads')
	plt.title('Quality Score Distribution over all Sequences')
	return plt

def __main__():
	all_ids = [] #! saved as object code, not as string! saved in item.name
	ascii = ("!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
	arguments = Parser()
	phred = arguments.phred
	count = 0

	with open(arguments.Dateipfad) as inhalt: #open user-provided file path
		line_pack = []
		test_cnt = 0	#!! KÜRZEN FÜR TESTLÄUFE
		for lines in inhalt:
			if line_pack == []:
				phred = Phred_Bestimmung(lines,ascii,phred)
			line_pack.append(lines.rstrip())
			test_cnt += 1	#!! KÜRZEN FÜR TESTLÄUFE
			if len(line_pack) == 4:
				id = lines[0]
				all_ids.append(Read(line_pack,id,ascii,phred,all_ids))
				line_pack = []
			if test_cnt > 10000:	#!! KÜRZEN FÜR TESTLÄUFE
				break

	fig = plt.figure()	#create canvas for graphs
	graph_gc(all_ids,fig)	#create first subplot, for gc-percentages
	graph_scores(all_ids,fig)	#fill in second subplot, for scores
	plt.show()
	# if arguments.save:
	# 	plt.savefig(str(arguments.save))

if __name__ == "__main__":
	__main__()
