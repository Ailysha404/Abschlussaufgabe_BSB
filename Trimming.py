scores=[176,299,398,409,5555,5,18,2,22,13,16,99,156,87,77] #Liste von Belinda
def trimming(scores):
    counter=0
    summe=0
    item_index=0
    kmer_neu=[]

    for kmer in scores[::-1]:
        counter+=1
        summe+=kmer
        item_index+=1
        if counter==3:
            kmer_mean=summe/3
            counter=0
            summe=0 #Summe wieder zurücksetzen!
            if kmer_mean <25:
                kmer_neu = scores[:len(scores)-item_index]
                return kmer_neu
            
print(trimming(scores))


        
    
    
    