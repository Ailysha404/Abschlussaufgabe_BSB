#Graphen

#GC Gehalt
import matplotlib.pyplot as plt
X=[0,1,2,3,4,5,6,7,8,9,10]#"Mean GC content(%)"
Y=[2**n for n in X]#"Count"
farbe=['blue']
plt.plot(X, Y, linewidth=2.0,alpha=0.7,label='Qualität der Basen')
plt.xlabel('Mean GC content (%)')
plt.ylabel('Count')
plt.title('GC distribution over all sequences')
plt.legend()
plt.show()

