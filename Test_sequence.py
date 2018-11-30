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

# def isGene3(seq): # return True if seq has start and stop codon
#     for i in range(len(seq)):
#         if isCodonStart(seq,i)==True:
#             for j in range(i,len(seq),3):
#                 if isCodonStop(seq,j)==True:
#                     print ("True")
#     return False

def isGene3(seq):
    for i in range(len(seq)-2):
        if isCodonStart(seq,i) == True:
            for j in range(i,len(seq)-2,3):
                if isCodonStop(seq,j) == True:
                    print("===================== Frame :",i%3,"===================== \n")
                    print("Length:",j-i,"pb")
                    print("Codon Start : "+oneWord(seq,i,3)+", Position : ",i,"; Codon Stop : "+oneWord(seq,j,3)+" position : ",j)
                    break
    return False

def openFasta(file):
    with open(file) as fasta:
        data = fasta.read()
        data = str.replace(data,"\n","")
    return data

if __name__=='__main__':
    data = openFasta("sequence.fasta")
    isGene3(data)
