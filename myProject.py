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
    '''
    La fonction créé un dictionnaire associant un codon en clé à sa traduction acide aminé en valeur

        Return :
                dictionnaire de la table génétique de Mycoplasma Genitallium
    '''
    dico={'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W','CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R','ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
    return dico

def trad(seq,cadre):
    '''
    La fonction prend une séquence et la traduit en acide aminé

        argument:
                seq:séquence Fasta utilisée
                cadre:cadre de lecture sur lequel traduire la séquence (1,2 ou 3)
        Return:
                la séquence traduite en format .fasta
    '''
    seqprot = ""
    for i in range(cadre,len(seq)-2,3):
        codon = ""
        prot = ""
        codon += seq[i]+seq[i+1]+seq[i+2]
        prot += dico_table[codon]
        seqprot += prot
    return seqprot

def openFasta(file):
    '''
    La fonction à pour objectif de récuperer le fichier Fasta présent dans le repertoire
        courant, elle le convertit en string, en enlevant les sauts de lignes et le stocke dans une variable

            argument:
                    file: nom du fichier source

                    return:
                    La fonction retourne data, c'est la variable qui contient la séquence Fasta standardisée pour le programme.
    '''
    with open(file) as fasta:
        data = fasta.readlines()
        data.pop(0)
        datastr = ""
        for i in data:
            datastr+=i
        datastr = str.replace(datastr,"\n","")
    return datastr

def writeFasta(data,FICHIER):
    '''
    Cette fonction n'as pas d'utilité à proprement parler dans l'éxécution du programme
        elle sert à écrire dans un fichier txt, les séquences obtenues à la fin de l'éxécution de four-lectures

            argument:
                    data:le dictionnaire où les valeurs sont placées et que vous voulez enregistrer
    '''
    fic=open(FICHIER, "w")
    for letter in data :
        fic.write(letter)
    fic.close()

def writeCSV(dictionary,filename,separator):
    '''
    Cette fonction enregistre certaines informations sur les ORFs contenues dans le 'dictionary' dans un fichier csv

            argument:
                    dictionary:le dictionnaire dont vous voulez enregistrer des données
                    filename: le nom du fichier .csv dans lequel les informations seront enregistrées
                    separator: le type de séparator utilisé pour créé votre fichier .csv (par défaut ",")
            return:
                    Rien
    '''
    listeCSV = []
    listeCSV.append("'cadre'"+separator+"'id'"+separator+"'start'"+separator+"'stop'"+"\n")
    for cadre in dictionary.keys():
        for orf in dictionary[cadre].keys():
            id, start, stop = orf, dictionary[cadre][orf]['Start'], dictionary[cadre][orf]['Stop']
            listeCSV.append(str(cadre)+separator+str(id)+separator+str(start)+separator+str(stop)+"\n")
    fic = open(filename,"w")
    for ligne in range(len(listeCSV)):
        fic.write(listeCSV[ligne])
    fic.close()

def readCSV(filename,separator):
    '''
    Cette fonction lit les informations contenues dans un fichier csv et les place dans un dictionnaire.

            argument:
                    filename: le nom du fichier .csv dans lequel les informations seront enregistrées
                    separator: le type de séparator utilisé pour créé votre fichier .csv (par défaut ",")
            return:
                    Le dictionnaire contenant les informations nouvellement récupérées.
    '''
    csvliste = []
    dicocsv = {}
    dicocsv[1] = {}
    dicocsv[2] = {}
    dicocsv[3] = {}
    orftmp=[]
    with open(filename) as csv:
        data = csv.readlines()
        for ligne in range(1,len(data)):
            csvliste.append(data[ligne])
            orftmp = csvliste[ligne-1].split(separator)
            orftmp[3] = orftmp[3][0:-1]
            taille = int(orftmp[3]) - int(orftmp[2])
            dicocsv[int(orftmp[0])][int(orftmp[1])] = {'Start':int(orftmp[2]),'Stop':int(orftmp[3]),'Taille (pb)':taille,'Seq_Nucleo':'NON RENSEIGNEE','Seq_proteo':'NON RENSEIGNEE'}
    return dicocsv

def fasta_csv_link(dicORF,seq):
    '''
    '''
    for cadre in dicORF.keys():
        for gene in dicORF[cadre]:
            seq_nucleot = seq[dicORF[cadre][gene]['Start']:dicORF[cadre][gene]['Stop']]
            dicORF[cadre][gene]['Seq_Nucleo'] = seq_nucleot
            dicORF[cadre][gene]['Seq_proteo'] = trad(seq_nucleot,0)
    return dicORF

def comp_reverse(seq):
    '''
    Cette fonction sert a inversé l'ordre de la séquence complémentaire.

        argument:
                seq: séquence Fasta utilisée

        return:
                La fonction retourne la séquence reverse complémentaire
    '''
    a = ''
    seq_comp = ''
    for i in seq:
        if i == 'A':
            a = 'T'
        elif i == 'T':
            a = 'A'
        elif i == 'G':
            a = 'C'
        elif i == 'C':
            a = 'G'
        seq_comp += a
    seq_comp_inv = seq_comp[::-1]
    return seq_comp_inv

def getLengths(orflist,display):
    '''
    La fonction parcourt le dictionnaire et renvoit la taille des ORFs dans une liste pour chaque cadre de lecture.

        argument:
                    orflist:dictionnaire dans lequel se trouve les informations

        Return:
                    Une liste contenant les longueurs des ORFs
    '''
    listelongueur=[]
    for cadre in orflist.keys():
        listelongueur.append([])
        for orf in orflist[cadre].keys():
            listelongueur[cadre-1].append(orflist[cadre][orf]['Taille (pb)'])
            if display == 0 :
                print("Pour le cadre : "+str(cadre)+" l'orf numéro : "+str(orf)+" a une taille de : "+str(orflist[cadre][orf]['Taille (pb)'])+"pb")
    return listelongueur

def getLongestORF(orflist):
    '''
    La fonction parcourt le dictionnaire et retourne les ORFs pour lesquels la taille est la plus importante

            argument:
                    orflist: dictionnaire dans lequel se trouve les informations

            Return:
                    Rien
    '''
    for cadre in orflist.keys():
        taillemax = 0
        for gene in orflist[cadre]:
            taille = orflist[cadre][gene]['Taille (pb)']
            if taille > taillemax :
                taillemax = taille
                cadremax = cadre
                genemax = gene
        print ("L'ORF ", genemax, " du cadre ", cadremax," avec une taille de : ", taillemax,"pb")

def getTopLongestORF(orflist,value):
    '''
    La fonction récupère la liste des longueurs des ORFs, les trie, ne garde qu'un pourcentage des tailles les plus importantes puis
        parcourt le dictionnaire contenant tous les ORFs et leurs caractéristiques et récupère les ORFs dont les tailles se trouvent dans
        la liste prédédemment récupérée. La fonction retourne enfin la liste de ces ORFs

            argument:
                    orflist:dictionnaire dans lequel se trouve les informations
                    value:pourcentage des ORFs les plus longs que l'on souhaite afficher

            Return:
                    Rien
    '''
    longestLengths = []
    lengths = getLengths(orflist,1)
    for i in range (3):
        lengths[i].sort()
        nombreORFcadre = len(lengths[i])
        if nombreORFcadre > 0:
            compteur = int(value*nombreORFcadre/100)
            if compteur == 0:
                compteur = 1
            for j in range(1,compteur+1):
                longestLengths.append(lengths[i][-j])
    for cadre in orflist.keys():
        for orf in orflist[cadre].keys():
            if orflist[cadre][orf]['Taille (pb)'] in longestLengths:
                print("Pour le cadre : "+str(cadre)+" l'orf numéro : "+str(orf)+" a une taille de : "+str(orflist[cadre][orf]['Taille (pb)'])+"pb")

def oneWord(seq,pos,wlen):
    '''
    La fonction sert à copier une séquence dans une liste à partir d'une position start et sur une longueur donnée

        argument:
                    wlen: longueur de la liste voulue
                    pos: position voulue de départ
                    seq:correspond a la séquence Fasta récupérer

        Return:     la fonction retourne une liste appelée séquence contenant
                    la séquence nucléotidique débutant à la position strat souhaitée et d'une longueur wlen
    '''
    a=0
    sequence=''
    for i in range(len(seq)):
        if i >= pos and a<wlen :
            a += 1
            sequence += str(seq[i])
    return sequence

def isCodonStart(seq,pos):
    '''
    La fonction sert à identifier le codon start de la séquence, en parcourant la séquence de 1 en 1 et stock dans la variable
        codon une suite de 3 nucléotides.

        argument:
                    seq:séquence Fasta
                    pos: position de départ du codon de 3 nucléotides à tester

        return:
                    La fonction retourne True si en parcourant la séquence,
                    elle trouve un des codons start de la liste, sinon elle retourne False
    '''
    codon=oneWord(seq,pos,3)
    if codon in ('ATG','TTA','TTG','CTG','ATT','ATC','ATA','GTG'):
        return True
    else :
        return False

def isCodonStop(seq,pos):
    '''
    la fonction sert à identifier le codon stop de la séquence, en parcourant la séquence de 1 en 1 et stock
            dans la variable codon une suite de 3 nucléotides

        argument:
                    seq: Séquence fasta
                    pos:position de départ du codon de 3 nucléotides à tester

        return:
                    La fonction retourne True, si en parcourant la séquence
                    elle trouve un des deux codons Stop, sinon elle retourne False
    '''
    codon=oneWord(seq,pos,3)
    if codon == 'TAA' or codon == 'TAG':
        return True
    else :
        return False

def isGene3(seq,version,threshold,fourchette,fork_basse,fork_haute,dicORF):
    '''
    La fonction parcourt la séquence sur sa longueur -2. Elle cherche alors en appellant la fonction isCodonStart
        si il y a présence d'un codon dans un premier cadre de lecture en parcourant de 3 en 3 la séquence. Si elle identifie un codon alors
        elle va alors appeler isCodonStop afin d'identifier un codon stop dans l'ORF associé au codon start identifié, elle s'arrète au premier dès qu'elle en trouve un. Si la distance du codon stop avec
        le codon start est inférieure au threshold, le programme repart du codon strat à la recherche d'un meilleur duo codon start/codon stop.
        Une fois un ORF déterminé, la fonction va print la version de la séquence étudiée (5' 3' ou complémentaire inversé),
        le cadre de lecture (0,1,2), la longueur j+2-i égale au début du codon start et à la fin du codon stop. De plus, elle print la suite nucléotidique du codon start et du codon stop
        et leur position. A la fin d'un cadre de lecture, elle reprend du début pour faire le deuxième cadre de lecture puis le troisième.

        arguement:
                    seq: séquence Fasta étudiée
                    version: la version peut être 5'3' ou complementaire inversé
                    threshold: cela correspond à une longueur d'ORF minimal afin de limiter les faux positifs
                    fourchette:
                                True: la fonction analyse la séquence limitée par une fouchette supérieure et inférieure
                                False: la fonction analyse toute la séquence
                    fork_basse: limite inférieure
                    fork_haute: limite supérieure

        return:
                    La fonction retourne un dictionnaire dans lequel toutes les informations des ORFs trouvées est inscrit
    '''
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

def menu(choix1,dico_setup,fichier_setup,dicoForward,dicoBackward):
    while choix1 != "CSV" and choix1 != "FASTA":
        print("Que voulez vous faire ?\n               Tapez 'CSV' pour lire les résultats d'une précédente étude.\n               Tapez 'FASTA' pour importer une séquence dans le programme au format fasta.\n")
        while True :
            try :
                choix1 = str(input())
                break
            except :
                print("truc")
        if choix1 == "CSV":
            print("Quel fichier voulez-vous lire ?")
            try :
                FICHIER = str(input())
            except ValueError:
                FICHIER = "truc.csv"
            dicoForward = readCSV(FICHIER,",")
            fichier_setup = True
            dico_setup = True
            print("Voulez vous ajouter un fichier fasta lié à votre fichier csv ? (o/n)")
            choix4 = str(input())
            if choix4 == "o":
                print("Quel fichier voulez-vous lire ?")
                try :
                    FICHIER = str(input())
                except ValueError :
                    FICHIER = "sequence.fasta"
                data = openFasta(FICHIER)
                data_inv_comp = comp_reverse(data)
                fichier_setup = True
                dicoForward = fasta_csv_link(dicoForward,data)
                # dicoBackward=fasta_csv_link(dicoBackward,data_inv_comp)
        elif choix1=="FASTA":
            print("Quel fichier voulez-vous lire ?")
            try :
                FICHIER = str(input())
            except ValueError :
                FICHIER = "sequence.fasta"
            data = openFasta(FICHIER)
            data_inv_comp = comp_reverse(data)
            fichier_setup = True
    print("Que voulez vous faire avec ces fichiers ?\n           Tapez 'AFFICHAGE' pour afficher la liste des orfs de votre fichier.\n           Tapez 'LONGEST' pour afficher les orfs les plus longs pour chaque cadre de lecture.\n           Tapez 'LONG' pour afficher une portion des orfs les plus longs.\n           Tapez 'WRITE' pour enregistrer vos résultats dans un fichier.\n           Tapez 'FIND ORF' pour trouver les orfs de votre séquence.\n           Tapez 'EXIT' pour quitter le programme !")
    choix2=str(input())
    if choix2 == 'AFFICHAGE':
        if dico_setup == True :
            print("Pour la séquence Forward : \n")
            getLengths(dicoForward,0)
            print("Pour la séquence Backward : \n")
            getLengths(dicoBackward,0)
        elif fichier_setup == True :
            print("Vous devez d'abord trouver les ORFs de votre fichier !")
            choix2 = 'FIND ORF'
    elif choix2 == 'LONGEST':
        if dico_setup == True :
            print("Voici la liste des ORFs les plus longs pour chaque cadre de lecture :")
            print("       Sens Forward")
            getLongestORF(dicoForward)
            print("       Sens Backward")
            getLongestORF(dicoBackward)
        elif fichier_setup == True :
            print("Vous devez d'abord trouver les ORFs de votre fichier !")
            choix2 = 'FIND ORF'
    elif choix2 == 'LONG':
        if dico_setup == True :
            print("Quelle proportion des ORFs les plus longs souhaitez-vous (en %)? (De base 50%)")
            try :
                pourcent = int(input())
            except :
                pourcent = 50
            print("Les ",pourcent,"% ORF les plus grands")
            print("Pour la séquence Forward : \n")
            getTopLongestORF(dicoForward,pourcent)
            print("Pour la séquence Backward : \n")
            getTopLongestORF(dicoBackward,pourcent)
        elif fichier_setup == True:
            print("Vous devez d'abord trouver les ORFs de votre fichier !")
            choix2='FIND ORF'
    elif choix2 == 'WRITE':
        if dico_setup == True :
            print("Quel est le nom du fichier sur lequel vous souhaitez enregistrer les données ?")
            try :
                FICHIER2 = str(input())
            except :
                FICHIER2 = "truc2.csv"
            if ".csv" in FICHIER2 :
                FICHIER2 = FICHIER2[0:-4]
            FICHIER2 = FICHIER2+"_forward.csv"
            FICHIER3 = FICHIER2[0:-12]+"_backward.csv"
            print("Les fichiers seront : ", FICHIER2," et ", FICHIER3)
            writeCSV(dicoForward,FICHIER2,",")
            writeCSV(dicoBackward,FICHIER3,",")
        elif fichier_setup == True :
            print("Vous devez d'abord trouver les ORFs de votre fichier !")
            choix2 = 'FIND ORF'
    if choix2 == 'FIND ORF':
        if dico_setup == True :
            print("Voulez vous changer de fichier à étudier ? (o/n)")
            while True :
                try :
                    choix3 = str(input())
                    if choix3 == "o" or choix3 == "n" :
                        break
                except :
                    break
            if choix3 == "o":
                choix1 = ""
        else :
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
            dicoForward = isGene3(data,"Forward",threshold,fourchette,fork_basse,fork_haute,dicoForward)
            dicoBackward = isGene3(data_inv_comp,"Backward",threshold,fourchette,fork_basse,fork_haute,dicoBackward)
            dico_setup=True
    elif choix2=='EXIT':
        return False,choix1,dico_setup,fichier_setup,dicoForward,dicoBackward
    return True,choix1,dico_setup,fichier_setup,dicoForward,dicoBackward


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
    choix1 = ''
    dico_setup = False
    fichier_setup = False
    go = True
    dicORF = {}
    dicORF_inv_comp = {}
    #On récupère un seuil pour la taille minimale des gènes
    while go == True :
        go,choix1,dico_setup,fichier_setup,dicORF,dicORF_inv_comp = menu(choix1,dico_setup,fichier_setup,dicORF,dicORF_inv_comp)
    ## On ouvre le fichier contenant le gène d'intérêt
    # data = openFasta("sequence.fasta")
    # writeFasta(trad(data,0),"fasta_prot.txt")
    # writeFasta(trad(data,1),"fasta_prot1.txt")
    # writeFasta(trad(data,2),"fasta_prot2.txt")
    # writeFasta(trad(data_inv_comp,0),"fasta_prot2.txt")
    # writeFasta(trad(data_inv_comp,1),"fasta_prot1.txt")
    # writeFasta(trad(data_inv_comp,2),"fasta_prot2.txt")
    # data_inv,data_comp,data_inv_comp=comp_reverse(data)
    # writeFasta(data_inv_comp,"fasta_inv_comp.txt")
    # writeFasta(data_inv,"fasta_inv.txt")
    # writeFasta(data_comp,"fasta_comp.txt")
    # dicoForward = isGene3(data,"Forward",threshold,fourchette,fork_basse,fork_haute)
    # print("---- Données sur le brin principal ----")
    # print(getLengths(dicoForward))
    # print(getLongestORF(dicoForward)) #A REVOIR
    # print(getTopLongestORF(dicoForward,50))
    # writeCSV(dicoForward,"ORFsearch.csv",";")
    # dicoBackward = isGene3(data_inv_comp,"Complémenaire Inverse",threshold,fourchette,fork_basse,fork_haute)
    # print("---- Données en complémentaire inverse ----")
    # print(getLengths(dicoBackward))
    # print(getLongestORF(dicoBackward))
    # print(readCSV("prout.csv",";"))
