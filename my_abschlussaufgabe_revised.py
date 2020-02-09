import sys
import argparse
import matplotlib.pyplot as plt

class Read: #class for all reads
	def __init__(self,content,id,phred,ascii,idlist):
		self.length = len(content[1])
		self.gc = (content[1].count("G")+content[1].count("C"))/self.length #in decimal - between 0 and 1
		self.name = content[0].split()[0] #how to get specific substring from line 1? #saves own id as name
		self.sequenz = content[1] #sequence-string
		self.quality = self.qualitaet(content,phred,ascii) #list of quality-scores

	def qualitaet(self,content,ascii,phred):
		values = {}
		#values[content[0]] = content[3] #keys ids, values score

		if phred == "33": #check for phred format
			alphabet_33 = {char:value-33 for value,char in enumerate(ascii,33) if value < 74} #make dictionary w ascii-character + phred score

			for score_string in alphabet_33: #id read
				score_list = [] #list with all scores from the read

				for letter in alphabet_33[score_string].rstrip(): #individual letter in score line
					score_list.append(alphabet_33[str(letter)]) #add score to list
					values[score_string.rstrip()] = score_list  #dictionary with id string as key, list with scores as value
					#values.append(alphabet_33[str(letter)])


		elif phred == "64":
			alphabet_64 = {char:value-64 for value,char in enumerate(ascii,33) if value > 63 and value < 105} #make dictionary w ascii-character + phred score

			for score_string in alphabet_64:
				score_list = []

				for letter in content[3]:
					#values.append(alphabet_64[str(letter)])
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
	input_parser.add_argument('-p','--phred', dest='opt', help="Optionale manuelle Angabe des Phred-Formats, 33 oder 64")
	input_parser.add_argument('-t','--trim', dest='trim_val', help="Optionale manuelle Angabe des Trimming-Scores",default=25)
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

def graph_gc(dict_reads):
	X=[] 	# GC content(%)
	#Y=[] 	# Count

	for value in dict_reads:
		X.append(value.gc*100)

	plt.hist(X,histtype="step",color="orange")
	plt.xlabel('GC content (%)')
	plt.ylabel('Count')
	plt.title('GC distribution over all sequences')
	plt.legend()
	plt.show()

def __main__():
	all_ids = [] #! saved as object code, not as string! saved in item.name
	ascii = ("!\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
	arguments = Parser()
	phred = arguments.opt
	count = 0

	with open(arguments.Dateipfad) as inhalt: #open user-provided file path
		line_pack = []
		#test_cnt = 0	#!! KÜRZEN FÜR TESTLÄUFE
		for lines in inhalt:
			count += 1
			print(count)
			if line_pack == []:
				phred = Phred_Bestimmung(lines,ascii,phred)
			line_pack.append(lines.rstrip())
			#test_cnt += 1
			if len(line_pack) == 4:
				id = lines[0]
				all_ids.append(Read(line_pack,id,ascii,phred,all_ids))
				line_pack = []
			#if test_cnt > 500:	#!! KÜRZEN FÜR TESTLÄUFE
			#	break
	count = 0

	for item in all_ids:
		count += 1
		print(count)
#		print(item.name) #Test length trimming - include buffer at the start?
#		print(len(item.quality[item]))
		item.quality = trimming(item.quality[item])
		item.sequenz = item.sequenz[0:len(item.quality)]
#		print(len(item.quality))
	graph_gc(all_ids)

if __name__ == "__main__":
	__main__()

# X=[0,1,2,3,4,5,6,7,8,9,10]#"Mean GC content(%)"
# Y=[2**n for n in X]#"Count"
# farbe=['blue']
# plt.plot(X, Y, linewidth=2.0,alpha=0.7,label='Qualität der Basen')
# plt.xlabel('Mean GC content (%)')
# plt.ylabel('Count')
# plt.title('GC distribution over all sequences')
# plt.legend()
# plt.show()
