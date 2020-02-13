import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import csv
from tabulate import tabulate

# '''Tabelle speichern in csv'''
test_table = {("a", 500, 500, 30), 
               ("b", 30, 40, 2),
               ("c", 40, 40, 5)}
print(tabulate((test_table), headers="keys"))
f = open("table.csv","w")
with f:
    writer = csv.writer(f,delimiter=';')
    for row in test_table:
        writer.writerow(row)
f.close()
# #https://docs.python.org/3/library/csv.html
# 
# '''csv.reader(csv-datei,Dialekt='exel',**formatparameter)'''
# # mit newline='' öfnnen wenn file object
# # Gibt ein Leserobjekt zurück, das über Zeilen in der angegebenen csv-Datei iteriert
# # Beliebiges Objekt mit __next__()-methode, Dateiobjekte und Listenobjekte sind beide geeignet
# '''Jede aus der csv-Datei gelesene Zeile wird als Liste von Zeichenketten zurückgegeben.
#     Es wird keine automatische Datentypkonvertierung durchgeführt, es sei denn, die
#     Formatierungsoption QUOTE_NONNUMERIC wird angegeben (in diesem Fall werden nicht in
#     Anführungszeichen gesetzte Felder in Fließkommazahlen umgewandelt).'''
# 
# # with open('eggs.csv', newline='') as csvfile:
# #     spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
# #     for row in spamreader:
# #         print(', '.join(row))
#         
# '''csv.writer(csv_datei,dialekt='excel',**formatparameter)'''
# # kann ein beliebiges Objekt mit einer write()-Methode
# # mit newline='' öffnen wenn file object
# '''wert None als leere Zeichenkette alle anderen Nicht-zeichenkettendaten als str()
#     vor dem schreiben mit einer Zeichenkette versehen'''
# '''Gibt ein Writer-Objekt zurück, das für die Konvertierung der Daten des Benutzers
#     in abgegrenzte Zeichenfolgen auf dem gegebenen dateiähnlichen Objekt sorgt'''
# 
# # 
# # with open('eggs.csv', 'w', newline='') as csvfile:
# #     spamwriter = csv.writer('eggs.csv', delimiter=' ',quotechar='|', quoting=csv.QUOTE_MINIMAL)
# #     spamwriter.writerow(['Spam'] * 5 + ['Baked Beans'])
# #     spamwriter.writerow(['Spam', 'Lovely Spam', 'Wonderful Spam'])
# 

 
     
'''Graphen'''
# '''Seaborn-Countplot'''
# # # seaborn countplot:
# # #     sns.countplot(x=readlänge)->Häufigkeit der Readlängen vorkommen
# # #     sns.distplot(readlänge,bins=50)
# # #     plt.show()
#         '''qual_plot=sns.distplot(y_mean,kdr=True,bins=50)#bennenung von graphen möglich,
#         für spätere direktausgabe über main'''

'''pandas.Dataframe-Boxplot'''#graph_gc(dict_reads, fig)
def gc_boxplot():

    df = pd.DataFrame((),columns=['','',''])

    boxplot = df.boxplot(column=['','',''])

    plt.show()

'''Lineplot-pandas'''
def basen_lineplot():
    df = pd.DataFrame({

       'C': [20, 18, 489, 675, 1776],

       'A': [4, 25, 281, 600, 1900],
       
       'T':[],
       
       'G':[]

       }, index=[Postition],label="Base Quality per Position")

    lines = df.plot.line()
    plt.show()

'''Countplot'''#Häufigkeit der Readlängen vorkommen
df = pd.melt(testdf)
sns.countplot(data=df.loc[df['value']!="NO"], x='variable', hue='value')
#->sns.countplot(x=readlänge)
plt.show()

'''Histogram mit pandas''' #graph_basequality(dict_reads, fig, dataframe)
plt.hist(df.x, alpha=.3)
sns.rugplot(df.x);
plt.show()

'''distplot mit pandas'''#raph_scores(dict_reads, fig)
sns.distplot(df.x)
plt.show()

#https://chrisalbon.com/python/data_wrangling/pandas_with_seaborn/