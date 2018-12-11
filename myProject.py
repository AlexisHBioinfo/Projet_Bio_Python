"""
22/11/2018
ANDRE Charlotte
CLAUDE Elsa
HUCTEAU Alexis
LAPORTE Antoine
BioPython
"""
from statistics import mean

def table_genetic():
    dico={'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W','CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R','ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
    return dico

def trad(seq,cadre):
    seqprot = ""
    for i in range(cadre,len(seq)-2,3):
        codon = ""
        prot = ""
        codon+=seq[i]+seq[i+1]+seq[i+2]
        prot+=dico_table[codon]
        seqprot+=prot
    return seqprot

def openFasta(file):
    '''La fonction à pour objectif de récuperer le fichier Fasta présent dans le repertoire
    courant, elle le convertit en string, en enlevant les sauts de lignes et le stocke dans une variable

    argument: file: nom du fichier source

    return: La fonction retourne data, c'est la variable qui contient la séquence Fasta standardisée pour le programme.
    '''
    with open(file) as fasta:
        data = fasta.readlines()
        data.pop(0)
        datastr=""
        for i in data:
            datastr+=i
        datastr = str.replace(datastr,"\n","")
    return datastr

def writeFasta(data,FICHIER):
    '''Cette fonction n'as pas d'utilité à proprement parler dans l'éxécution du programme
    elle sert à écrire dans un fichier txt, les séquences obtenues à la fin de l'éxécution de four-lectures'''
    fic=open(FICHIER, "w")
    for letter in data :
        fic.write(letter)
    fic.close()

def writeCSV(dictionary,filename,separator):
    listeCSV=[]
    for cadre in dictionary.keys():
        for orf in dictionary[cadre].keys():
            id, start, stop = orf, dictionary[cadre][orf]['Start'], dictionary[cadre][orf]['Stop']
            listeCSV.append("'cadre'"+str(cadre)+separator+"'id' :"+str(id)+separator+"'start'"+str(start)+separator+"'stop'"+str(stop)+"\n")
    fic = open(filename,"w")
    for ligne in range(len(listeCSV)):
        fic.write(listeCSV[ligne])
    fic.close()

def readCSV(filename,separator):
    csvliste=[]
    with open(filename) as csv:
        data = csv.readlines()
        for ligne in range(len(data)):
            csvliste.append(data[ligne])
            print(csvliste[ligne])
    return csvliste

def four_lectures(seq):
    '''Cette fonction sert a inversé l'ordre de la séquence, obtenir la séquence complémentaire, et inversé l'ordre de la séquence complémentaire.

    argument: seq: séquence Fasta utilisée

    return: La fonction return les 4 séquences obtenues
    '''
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

def getLengths(orflist):
    listelongueur=[]
    for cadre in orflist.keys():
        listelongueur.append([])
        for orf in orflist[cadre].keys():
            listelongueur[cadre-1].append(orflist[cadre][orf]['Taille (pb)'])
            print("Pour le cadre : "+str(cadre)+" l'orf numéro : "+str(orf)+" a une taille de : "+str(orflist[cadre][orf]['Taille (pb)'])+"pb")
    return listelongueur

def getLongestORF(orflist):
    maxval = []
    for i in range(1,4):
        print(i)
        try:
            maxval.append(max(orflist[i].items(),key=lambda x:x[1]['Taille (pb)']))
            print("Pour le cadre",i,"l'ORF",maxval[0],"a la taille maximale de",maxval[1]['Taille (pb)'])
        except ValueError:
            print("Pas de données pour le cadre",i)
    return maxval

def getTopLongestORF(orflist,value):
    longestLengths = []
    ListeLongestORFLength =[]
    lengths = getLengths(orflist)
    for i in range (3):
        lengths[i].sort()
        compteur = int((value / 100) * len(lengths[i]))
        for j in range(1,compteur+1):
            longestLengths.append(lengths[i][-j])
    for cadre in orflist.keys():
        for orf in orflist[cadre].keys():
            if orflist[cadre][orf]['Taille (pb)'] in longestLengths:
                ListeLongestORFLength.append("Pour le cadre : "+str(cadre)+" l'orf numéro : "+str(orf)+" fait parti des "+str(value)+"% orf les plus grands, avec une taille de : "+str(orflist[cadre][orf]['Taille (pb)'])+"pb")
    return ListeLongestORFLength

#2.2.1 Statistiques basiques

def oneWord(seq,pos,wlen): #return a piece with a lenght = wlen of the string seq from the index start
    '''
    La fonction sert à copier une séquence dans une liste à partir d'une position start et sur une longueur donnée

        argument: wlen: longueur de la liste voulue
        pos: position voulue de départ
        seq:correspond a la séquence Fasta récupérer

        Return: la fonction retourne une liste appelée séquence contenant
        la séquence nucléotidique débutant à la position strat souhaitée et d'une longueur wlen
    '''
    a=0
    sequence=''
    for i in range(len(seq)):
        if i>=pos and a<wlen :
            a=a+1
            sequence=sequence+str(seq[i])
    return sequence

def isCodonStart(seq,pos): #return True if the string seq has start codon
    '''
    La fonction sert à identifier le codon start de la séquence, en parcourant la séquence de 1 en 1 et stock dans la variable
        codon une suite de 3 nucléotides.

        argument: seq:séquence Fasta
                pos: position de départ du codon de 3 nucléotides à tester

        return: La fonction retourne True si en parcourant la séquence,
        elle trouve un des codons start de la liste, sinon elle retourne False
    '''
    codon=oneWord(seq,pos,3)
    if codon in ('ATG','TTA','TTG','CTG','ATT','ATC','ATA','GTG'):
        return True
    else :
        return False

def isCodonStop(seq,pos): #return True if seq has a stop codon
    '''
    la fonction sert à identifier le codon stop de la séquence, en parcourant la séquence de 1 en 1 et stock
            dans la variable codon une suite de 3 nucléotides

        argument: seq: Séquence fasta
                pos:position de départ du codon de 3 nucléotides à tester

        return: La fonction retourne True, si en parcourant la séquence
        elle trouve un des deux codons Stop, sinon elle retourne False
    '''
    codon=oneWord(seq,pos,3)
    if codon=='TAA' or codon=='TAG':
        return True
    else :
        return False

def isGene3(seq,version,threshold,fourchette,fork_basse,fork_haute):
    '''
    La fonction parcours la séquence sur sa longueur -2. Elle cherche alors en appellant la fonction isCodonStart
        si il y a présence d'un codon sans cadre de lecture spécifique. Si elle identifie un codon alors
        elle parcourt de 3 en 3 la séquence en commencant à la position i+threshold, qui correspond a la position
        du codon start plus un seuil minimal d'ORF,par défaut cette valeur de threshold est de 90. Une fois cette distance atteinte,
        la fonction va alors appeler isCodonStop afin d'identifier un codon stop dans l'ORF associé au codon start identifié, elle s'arrète au premier dès qu'elle en trouve un.
        Une fois un ORF déterminé, la fonction va print la version de la séquence étudiée (5' 3' ou complémentaire inversé),
        le cadre de lecture (0,1,2), la longueur j+2-i égale au début du codon start et à la fin du codon stop. De plus, elle print la suite nucléotidique du codon start et du codon stop
        et leur position, sinon elle retourne False

        arguement: seq: séquence Fasta étudiée
                    version: la version peut être 5'3' ou complementaire inversé
                    threshold: cela correspond à une longueur d'ORF minimal afin de limiter les faux positifs

        return: La fonction ne retourne rien, mais elle print des informations
    '''
    dicORF = {}
    for k in range(1,4):
        dicORF[k] = {}
        listeGenes = []
        if fourchette == False:
            i = k-1
            lenSeq = len(seq)
        else:
            i = k-1+fork_basse
            lenSeq = fork_haute
        compteur = 0
        while i <= (lenSeq-2):
            if (i+1-k)%150 == 0:
                print("Le programme en est à la position :",i+1, "frame :",k)
            if isCodonStart(seq,i) == True:
                for j in range(i,lenSeq-2,3):
                    if isCodonStop(seq,j) == True and j+3-i >= threshold:
                        compteur+=1
                        print("=====================",version,": Frame :",k,"===================== ")
                        print("Length:",j+3-i,"pb")
                        print("Codon Start : "+oneWord(seq,i,3)+", Position : ",i+1,"; Codon Stop : "+oneWord(seq,j,3)+" position : ",j+3,"\n")
                        listeGenes.append("Codon Start : "+oneWord(seq,i,3)+", Position : "+str(i+1)+"; Codon Stop : "+oneWord(seq,j,3)+" position : "+str(j+3))
                        listeGenes.append("Sequence nucleotidique : "+oneWord(seq,i,j+3-i))
                        listeGenes.append("Sequence proteique : "+trad(oneWord(seq,i,j+3-i),0)+"\n")
                        dicORF[k][compteur]={'Start':i+1,'Stop':j+3,'Taille (pb)':j+3-i,'Seq_Nucleo':oneWord(seq,i,j+3-i),'Seq_proteo':trad(oneWord(seq,i,j+3-i),0)}
                        i=j
                        break
                    elif isCodonStop(seq,j) == True and j+3-i < threshold:
                        break
            i+=3
        fichierGene=""
        for i in listeGenes:
            fichierGene+=i+"\n"
        writeFasta(fichierGene,"cadre"+str(k)+".txt")
    return dicORF

def menu():
    choix1=''
    dico_setup=False
    fichier_setup=False
    while choix1!="CSV" and choix1!="FASTA":
        print("Que voulez vous faire ?\n       Tapez 'CSV' pour lire les résultats d'une précédente étude.\n       Tapez 'FASTA' pour importer une séquence dans le programme au format fasta.\n")
        while True :
            try :
                choix1=str(input())
                break
            except :
                print("truc")
    if choix1=="CSV":
        print("Quel fichier voulez-vous lire ?")
        try :
            FICHIER=str(input())
        except ValueError:
            FICHIER="truc.csv"
        readCSV(FICHIER,";")
        fichier_setup=True
        dico_setup=True
    elif choix1=="FASTA":
        print("Quel fichier voulez-vous lire ?")
        try :
            FICHIER=str(input())
        except ValueError :
            FICHIER="sequence.fasta"
        data = openFasta(FICHIER)
        fichier_setup=True
    print("Que voulez vous faire avec ces fichiers ?\n       Tapez 'AFFICHAGE' pour afficher la liste des orfs de votre fichier.\n       Tapez'LONGEST' pour afficher les orfs les plus longs pour chaque cadre de lecture.\n       Tapez 'LONG' pour afficher une portion des orfs les plus longs.\n       Tapez 'WRITE' pour enregistrer vos résultats dans un fichier.\n       Tapez 'FIND ORF' pour trouver les orfs de votre séquence.")
    choix2=str(input())
    if choix2=='AFFICHAGE':
        while True:
            if dico_setup==True :
                getLengths(dicoForward)
                break
            elif fichier_setup==True :
                print("Vous devez d'abord trouver les ORFs de votre fichier !")
                try:
                    threshold = int(input("Quelle valeur de threshold voulez-vous ? Entrez un nombre et appuyez sur ENTREE pour valider votre choix. La valeur par défaut est de 90pb : "))
                except ValueError:
                    threshold = 90
                    #On récupère une fourchette d'exécution du programme au cas où l'utilisateur souhaiterait cibler un endroit particulier
                try:
                    fourchette = str(input("Entrez une fourchette de lecture du génome (deux nombres séparés par un espace) et appuyez sur ENTREE, par défaut le génome entier est analysé : "))
                    fork_basse = (int(fourchette.split(" ")[0])//3)*3
                    fork_haute = (int(fourchette.split(" ")[1])//3)*3
                    if fork_basse > fork_haute: #On vérifie que l'utilisateur ne se soit pas trompé de sens dans l'entrée de la fourchette de lecture
                        fork_basse, fork_haute = fork_haute, fork_basse
                except ValueError:
                    fourchette = False
                    fork_basse = None
                    fork_haute = None
                dicoForward = isGene3(data,"Forward",threshold,fourchette,fork_basse,fork_haute)
                dico_setup=True
    elif choix2=='LONGEST':
        if dico_setup==True :
            print(getLongestORF(dicoForward))
        elif fichier_setup==True :
            choix2=='FIND ORF'
    elif choix2=='LONG':
        if dico_setup==True :
            print("Quelle proportion des ORFs les plus longs souhaitez-vous (en %)? (De base 50%)")
            try :
                porcent=str(input())
            except :
                porcent=50
            print(getTopLongestORF(dicoForward,porcent))
        elif fichier_setup==True:
            choix2=='FIND ORF'
    elif choix2=='WRITE':
        if dico_setup==True :
            print("Quel est le nom du fichier sur lequel vous souhaitez enregistrer les données ?")
            try :
                FICHIER2=str(input())
            except :
                FICHIER2="truc2.csv"
        elif fichier_setup==True :
            choix2=='FIND ORF'
    elif choix2=='FIND ORF':
        if dico_setup==True :
            print("Voulez vous changer de fichier à étudier ? (o/n)")
            while True :
                try :
                    choix3=str(input())
                    if choix3=="o" or choix3=="n" :
                        break
                except :
                    break
            if choix3=="o":
                choix1=""


######### WARNING ###########
#getGeneticCode(NCBI_ID) pour les biais de codon entre les espèces
#findORF(seq, threshold,codeTable) = isGene3 mais prend en compte les biais de Codon


#Ce serait bien de save les séquences inv, comp et comp inv dans un fasta ? DONE
# -> Ajout des codons init' alternatifs
# -> Suppression d'un codon stop (non présent dans le biais de codon)
# -> Ajout Threshold (minimum length : 90pb) "various thresholds (No threshold, 90bp, 210bp, 300bp, 420bp, for example.)"
#############################

if __name__=='__main__':
    dico_table = table_genetic()
    #On récupère un seuil pour la taille minimale des gènes
    menu()
    #On ouvre le fichier contenant le gène d'intérêt
    # data = openFasta("sequence.fasta")
    # writeFasta(trad(data,0),"fasta_prot.txt")
    # writeFasta(trad(data,1),"fasta_prot1.txt")
    # writeFasta(trad(data,2),"fasta_prot2.txt")
    # writeFasta(trad(data_inv_comp,0),"fasta_prot2.txt")
    # writeFasta(trad(data_inv_comp,1),"fasta_prot1.txt")
    # writeFasta(trad(data_inv_comp,2),"fasta_prot2.txt")
    # data_inv,data_comp,data_inv_comp=four_lectures(data)
    # writeFasta(data_inv_comp,"fasta_inv_comp.txt")
    # writeFasta(data_inv,"fasta_inv.txt")
    # writeFasta(data_comp,"fasta_comp.txt")
    #dicoForward = isGene3(data,"Forward",threshold,fourchette,fork_basse,fork_haute)
    # print("---- Données sur le brin principal ----")
    # print(getLengths(dicoForward))
    # print(getLongestORF(dicoForward)) #A REVOIR
    # print(getTopLongestORF(dicoForward,50))
    #writeCSV(dicoForward,"ORFsearch.csv",";")
    # dicoBackward = isGene3(data_inv_comp,"Complémenaire Inverse",threshold,fourchette,fork_basse,fork_haute)
    # print("---- Données en complémentaire inverse ----")
    # # print(getLengths(dicoBackward))
    # print(getLongestORF(dicoBackward))
    # print(readCSV("prout.csv",";"))
