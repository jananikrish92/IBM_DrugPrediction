import csv
import numpy as np
from calculate_psc import calculate_psc
def calPscIdx():
    with open('../BindingDB_All.tsv') as tsvfile:   
        reader=csv.reader(tsvfile, delimiter='\t')
        count=0
        greater_49=0
        less_49=0
        c=0 
        dict_prot={}
        for row in reader:
            
           # print(row)
            if(len(row)>49 and count!=1):
                greater_49+=1
                continue
                
            if(len(row)<49):
                less_49+=1
                continue
                
            if(row[41]!=''):
                #row[37]=
                dict_prot[row[41].split(',')[0]]=row[37]
          #  print row[41]
            count+=1

    psc_array = np.zeros([len(dict_prot),8420])
    j=0
    prot_id2idx = {}
    #foo=0
    for i in dict_prot.keys():
     #   foo+=1
      #  print dict_prot[i]
        psc_array[j]=calculate_psc(dict_prot[i])
        prot_id2idx[i] = j
        j+=1

    return psc_array,prot_id2idx
    #p.savetxt("foo.csv", g, delimiter=",")
  #  fd.write(g)
    #print g[0]
np.savetxt("psc.csv", psc_array.T, delimiter=",", header=','.join(dict_prot.keys()))