def calculation():
    '''Calculation of the mean,
    standard deviation and median for each read of the fastq file'''
    
    np.mean('''liste der einzelnen Reads''')
    np.median('''liste der einzelnen Reads''')
    np.std('''liste der einzelnen Reads''')
    
    #Speichern in einer Tabelle inklusive durchschnittlicher Readlänge zu dem Zeitpunkt
    pass #Tabelle
    return (||||<tablewidth="80%">'''Calculation_Data'''||
    ||<:> Mean ||<:>Median ||<:>Standard deviation ||<:>Timing||<:>Average read length||)

    '''Durchführen vor dem Filtern, nachdem Filtern(vor dem Trimmen)
    und nach dem Trimmen für GC_Gehalt und Verteilung der Readlängen
    (immer mit Durchschnittlicher Readlänge)'''
    
    return 

def graph_basequality(dict_reads,score_list,fig):
"""Shows the Quality of all Read positions"""
    fig.add_subplot(4,2,(1,4))
    X=[score_list]#x and y have to be the same lenght
    Y=[range(0,len(score_list))]
    
    plt.plot(Y,X,alpha=0.7,color="black")
    plt.xlabel('Posittion')
    plt.ylabel('Quality')
    plt.title('Quality Distribution over all Reads')
    
    return plt

    graph_basequality(all_ids,score_list,fig)