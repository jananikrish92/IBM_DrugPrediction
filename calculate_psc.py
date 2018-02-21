#a="MPPYISAFQAAYIGIEVLIALVSVPGNVLVIWAVKVNQALRDATFCFIVSLAVADVAVGALVIPLAILINIGPQTYFHTCLMVACPVLILTQSSILALLAIAVDRYLRVKIPLRYKTVVTQRRAAVAIAGCWILSLVVGLTPMFGWNNLSVVEQDWRANGSVGEPVIKCEFEKVISMEYMVYFNFFVWVLPPLLLMVLIYLEVFYLIRKQLNKKVSASSGDPQKYYGKELKIAKSLALILFLFALSWLPLHILNCITLFCPTCQKPSILIYIAIFLTHGNSAMNPIVYAFRIHKFRVTFLKIWNDHFRCQPKPPIDEDLPEEKAED"
#a="MDAMKRGLCCVLLLCGAVFVSPSQEIHARFRRGARSYQVICRDEKTQMIYQQHQSWLRPVLRSNRVEYCWCNSGRAQCHSVPVKSCSEPRCFNGGTCQQALYFSDFVCQCPEGFAGKCCEIDTRATCYEDQGISYRGTWSTAESGAECTNWNSSALAQKPYSGRRPDAIRLGLGNHNYCRNPDRDSKPWCYVFKAGKYSSEFCSTPACSEGNSDCYFGNGSAYRGTHSLTESGASCLPWNSMILIGKVYTAQNPSAQALGLGKHNYCRNPDGDAKPWCHVLKNRRLTWEYCDVPSCSTCGLRQYSQPQFRIKGGLFADIASHPWQAAIFAKHRRSPGERFLCGGILISSCWILSAAHCFQERFPPHHLTVILGRTYRVVPGEEEQKFEVEKYIVHKEFDDDTYDNDIALLQLKSDSSRCAQESSVVRTVCLPPADLQLPDWTECELSGYGKHEALSPFYSERLKEAHVRLYPSSRCTSQHLLNRTVTDNMLCAGDTRSGGPQANLHDACQGDSGGPLVCLNDGRMTLVGIISWGLGCGQKDVPGVYTKVTNYLDWIRDNMRP"
import numpy as np
dict_single={}
sampleseq="ARNDCEQGHILKMFPSTWYV"
for i in range(len(sampleseq)):
    dict_single[sampleseq[i]]=0


def calculate_psc(fasta_seq):
    
    fasta_seq=fasta_seq.upper()
   # dict_single={}
    for i in range(len(sampleseq)):
        dict_single[sampleseq[i]]=0
   # count=0
    for i in range(len(fasta_seq)):
        #fasta_seq[i]=fasta_seq[i]
        try:
            dict_single[fasta_seq[i]]+=1
        except:
          #  print foo
            empty_psc=np.zeros(8420)
            print "fail"
            return empty_psc
    #    count+=1

    for i in dict_single.keys():
        dict_single[i]=float(dict_single[i])/float(len(fasta_seq))
        
    dict_di={}
    for i in dict_single.keys():
        for j in dict_single.keys():
            dict_di[i+j]=0
    dict_tri={}
    for i in dict_single.keys():
        for j in dict_single.keys():
            for k in dict_single.keys():
                dict_tri[i+j+k]=0
    
    for i in range(len(fasta_seq)-1):
        dict_di[fasta_seq[i]+fasta_seq[i+1]]+=1
    
    for i in dict_di.keys():
        dict_di[i]=float(dict_di[i])/float(len(fasta_seq)-1)
        
    for i in range(len(fasta_seq)-2):
        dict_tri[fasta_seq[i]+fasta_seq[i+1]+fasta_seq[i+2]]+=1
    
    for i in dict_tri.keys():
        dict_tri[i]=float(dict_tri[i])/float(len(fasta_seq)-2)
    psc=np.zeros(len(dict_single)+len(dict_di)+len(dict_tri))
   # count=0
    psc[:len(dict_single)]=np.asarray(dict_single.values())
    psc[len(dict_single):len(dict_di)+len(dict_single)]=np.asarray(dict_di.values())
    psc[len(dict_di):len(dict_tri)+len(dict_di)]=np.asarray(dict_tri.values())
    #print psc.shape
    return psc