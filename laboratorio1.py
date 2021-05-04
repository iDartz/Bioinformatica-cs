from Bio import SeqIO
import matplotlib.pyplot as plt


def compara(seq1,seq2,t):
    auxt=0
    for i in range(len(seq1)):
        if seq1[i]==seq2[i]:
            auxt+=1
    
    if auxt>=t:
        return 1
    else:
        return 0

def Matrix(seq1,seq2,w,t):
    
    subseq1=[]
    subseq2=[]
    cadtem=""
    contador=0
    for i in seq1:
        if contador==w:
            subseq1.append(cadtem)
            cadtem=""
            contador=0
        else:
            contador+=1
            cadtem=cadtem+i
    
    cadtem=""
    contador=0
    for i in seq2:
        if contador==w:
            subseq2.append(cadtem)
            cadtem=""
            contador=0
        else:
            contador+=1
            cadtem=cadtem+i
    M = []
    
    for i in range(len(subseq1)):
        M.append([])
        for j in range(len(subseq2)):
            if compara(subseq1[i],subseq2[j],t)==1:
                M[i].append(1)
            else:
                M[i].append(0)
        
    return M

def plotMatrix(M):
    xs=[]
    ys=[]
    for i in range(len(M)):
        for j in range(len(M[i])):
            if M[i][j]==1:
                xs.append(i)
                ys.append(j)
    fig,ax=plt.subplots()
    ax.plot(xs,ys,'g.')
    plt.show()

def dotplot(seq1,seq2,w,t):
     
    
    M=Matrix(seq1,seq2,w,t)
    plotMatrix(M)



sequences = SeqIO.parse("P21333.fasta","fasta")
for record in sequences:
    data1=str(record.seq.upper())
    
sequences = SeqIO.parse("Q8BTM8.fasta","fasta")
for record in sequences:
    data2=str(record.seq.upper())
sequences = SeqIO.parse("Q8VHX6.fasta","fasta")
for record in sequences:
    data3=str(record.seq.upper())
    
sequences = SeqIO.parse("Q14315.fasta","fasta")
for record in sequences:
    data4=str(record.seq.upper())
threshold=23 #porcentaje
window=10

threshold=threshold*window/100

dotplot(data2,data1,window,threshold)
dotplot(data4,data3,window,threshold)
