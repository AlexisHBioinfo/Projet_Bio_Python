"""
22/11/2018
HUCTEAU
Alexis
BioPython Un peu plus loin avec Python (version déroulée)
"""

#2.2.1 Statistiques basiques


def isDNA(seq): #Return True if the string seq is a DNA sequence
    for i in seq :
        if not (i== "A" or i=="C" or i== "G" or i== "T" or i== "a" or i=="c" or i== "g" or i== "t") :
            return False
    return True

def countPro(seqProt) : #return the number of P in the string seqProt
    a=0
    for i in seqProt :
        if i == "P" or i == "p" :
            a=a+1
    return a

def countAll(seqProt,symbol): #return the number of symbol in the string seqProt
    a=0
    for i in seqProt :
        if i == symbol :
            a=a+1
    return a


def oneWord(seq,start,wlen): #return a piece with a lenght = wlen of the string seq from the index start
    a=0
    sequence=''
    for i in range(len(seq)):
        if i>=start and a<wlen :
            a=a+1
            sequence=sequence+str(seq[i])
    return sequence

def countWord(seq,word): #return the number of word in the string seq
    wlen=len(word)
    a=0
    for i in range(len(seq)):
        testlist=oneWord(seq,i,wlen)
        if testlist==word :
            a=a+1
    return a

def isCodonStart(seq,pos): #return True if the string seq has start codon
    codon=oneWord(seq,pos,3)
    if codon=='ATG':
        return True
    else :
        return False

def isCodonStop(seq,pos): #return True if seq has a stop codon
    codon=oneWord(seq,pos,3)
    if codon=='TAA' or codon=='TAG' or codon=='TGA' :
        return True
    else :
        return False

def isGene(seq): #return True if seq has start codon and stop codon modulo 3
    wlen=len(seq)
    for i in range(0,wlen,3):
        if isCodonStart(seq,i)==True :
            for j in range(i,wlen,3):
                if isCodonStop(seq,j)==True :
                    return True
    return False

def isGene3(seq): # return True if seq has start and stop codon
    for i in range(len(seq)):
        if isCodonStart(seq,i)==True:
            for j in range(i,len(seq),3):
                if isCodonStop(seq,j)==True:
                    print ("True")
    return False



"""
#########################isDNA FONCTION######################
seq='actgcgtaccttg'
flag = isDNA(seq) # False. This sequence contains a 'p'
print (flag)
#########################countPro###########################
seq='MVHLSAEEKEAVLGLWGKVNVDEVGGEALGRLLVVYPWTQ'
num = countPro(seq) # 1
print (num)
#########################countAll############################
seq='MVHLSAEEKEAVLGLWGKVNVDEVGGEALGRLLVVYPWTQ'
num = countAll(seq,'P') # 1
print (num)
num = countAll(seq,'V') # 7
print (num)
#########################oneWorld############################
seq='ABCDEFGHIJKLM'
s = oneWord(seq,4,5) # "EFGHI"
#########################countWord###########################
phrase='LESCHAUSSETTESDELARCHIDUCHESSESONTELLESSECHES'
num = countWord(phrase,'SSE') # 3
print (num)
#########################isCodonStart########################
seq='CTGATGTTCCATTACCAGTACAACAAACTATGATTCCATTACCAGTACA'
flag = isCodonStart(seq,3) # True
print(flag)
#########################isCodonStop#########################
seq='CTGATGTTCCATTACCAGTACAACAAACTATGATTCCATTACCAGTACA'
flag = isCodonStop(seq,30) # True
print(flag)
##########################isGene#############################
seq='ATGTGACTGATGTTCCATTACCAGTACAACAAACTATGATTCCATTACCAGTACA'
flag = isGene3(seq) # True
print(flag)
"""
