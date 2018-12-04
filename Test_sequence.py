"""
22/11/2018
ANDRE Charlotte
CLAUDE Elsa
HUCTEAU Alexis
LAPORTE Antoine
BioPython
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
    if codon in ('ATG','TTA','TTG','CTG','ATT','ATC','ATA','GTG'):
        return True
    else :
        return False

def isCodonStop(seq,pos): #return True if seq has a stop codon
    codon=oneWord(seq,pos,3)
    if codon=='TAA' or codon=='TAG':
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

def isGene3(seq,version,threshold=90):
    for i in range(len(seq)-2):
        if isCodonStart(seq,i) == True:
            for j in range(i+threshold,len(seq)-2,3):
                if isCodonStop(seq,j) == True:
                    #if j-i>=threshold :
                    print("=====================",version,": Frame :",i%3,"===================== ")
                    print("Length:",j-i,"pb")
                    print("Codon Start : "+oneWord(seq,i,3)+", Position : ",i+1,"; Codon Stop : "+oneWord(seq,j,3)+" position : ",j+1,"\n")
                    break
    return False

######### WARNING ###########
#getGeneticCode(NCBI_ID) pour les biais de codon entre les espèces
#findORF(seq, threshold,codeTable) = isGene3 mais prend en compte les biais de Codon


#Ce serait bien de save les séquences inv, comp et comp inv dans un fasta ? DONE
# -> Ajout des codons init' alternatifs
# -> Suppression d'un codon stop (non présent dans le biais de codon)
# -> Ajout Threshold (minimum length : 90pb) "various thresholds (No threshold, 90bp, 210bp, 300bp, 420bp, for example.)"
#############################

def openFasta(file):
    with open(file) as fasta:
        data = fasta.read()
        data = str.replace(data,"\n","")
    return data

def writeFasta(data,FICHIER):
    fic=open(FICHIER, "w")
    for letter in data :
        fic.write(letter)
    fic.close()

def four_lectures(seq):
    a=''
    seq_inv=seq[::-1]
    seq_comp=''
    for i in seq:
        if i =='A':
            a='T'
        elif i=='T':
            a='A'
        elif i=='G':
            a='C'
        elif i=='C':
            a='G'
        seq_comp += a
    seq_comp_inv=seq_comp[::-1]
    return seq_inv,seq_comp,seq_comp_inv

if __name__=='__main__':
    print ("Quelle valeur de threshold voulez-vous ? Entrez un nombre et tappez sur ENTREE pour valider votre choix. La valeur par défaut est de 90pb.")
    threshold=int(input())
    data = openFasta("sequence.fasta")
    data_inv,data_comp,data_inv_comp=four_lectures(data)
    # writeFasta(data_inv_comp,"fasta_inv_comp.txt")
    # writeFasta(data_inv,"fasta_inv.txt")
    # writeFasta(data_comp,"fasta_comp.txt")
    isGene3(data,"5'-3'",threshold)
    isGene3(data_inv,"3'-5'",threshold)
    isGene3(data_comp,"comp_5'-3'",threshold)
    isGene3(data_inv_comp,"comp_3'-5'",threshold)
