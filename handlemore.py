# using MODEST - a variant of IGLOSS with filtering - on one or more motifs
# motifs should be in files motif1.txt, motif2.txt, etc. 

proteome="AT.fasta"     # proteome in fasta format
numotifs=3              # number of motifs
scale=[10,10,10]        # list of scales, one for each motif
th=0                    # threshold 0=IGLOSS, usually 0.001

import os

allP=[]
allmotif=[]
allstart=[]
allstop=[]
allEvalue=[]
    
for k in range(1,numotifs+1):
    os.system("./MODEST "+ proteome + " motif"+str(k)+".txt " + str(scale[k-1]) + " 10 " + str(th) + " > " + " output.txt")

    f=open("output.txt","r")
    lines=f.readlines()
    f.close()
        
    for i in range(min(6,len(lines))):
        print(lines[i].rstrip())

    P=[]
    motif=[]
    start=[]
    stop=[]
    Evalue=[]
    ok=0
    for i in range(len(lines)):
        if "Evalue" in lines[i]:
            tempstart=int(lines[i].split()[1].split("-")[0])
            tempstop=int(lines[i].split()[1].split("-")[1])
            tempEvalue=float(lines[i].split("=")[1].split(")")[0])
        if ">" in lines[i]:
            tmp=lines[i].split(">")[1].split()[0]
            if tmp not in P:
                P.append(tmp)
                start.append(tempstart)
                stop.append(tempstop)
                Evalue.append(tempEvalue)
                ok=1
        if ok==1 and "Evalue" not in lines[i] and ">" not in lines[i]:
            motif.append(lines[i].rstrip())
            ok=0
            
    for i in range(len(P)):
        P[i]=P[i].split(".")[0]

    allP.append(P)
    allmotif.append(motif)
    allstart.append(start)
    allstop.append(stop)
    allEvalue.append(Evalue)


ind=[]
for i in range(len(allP[0])):
    ok=1
    indtemp=[]
    for k in range(1,numotifs):
        if allP[0][i] not in allP[k]:
            ok=0
    if ok==1:
        indtemp.append(i)
        for k in range(1,numotifs):
            indtemp.append(allP[k].index(allP[0][i]))
    if len(indtemp)>0:
        ind.append(indtemp)

os.system("rm output.txt")

f=open("MODESToutput.tsv","w")

f.write("SeqID\t")
for k in range(numotifs):
    f.write("motif"+str(k+1)+"\t"+"m"+str(k+1)+"_position_start\t"+"m"+str(k+1)+"_position_end\t"+"m"+str(k+1)+"_Evalue\t")
f.write("\n")
for i in range(len(ind)):
    f.write("%s\t" % allP[0][ind[i][0]])
    for k in range(numotifs):
        f.write("%s\t%d\t%d\t%f\t" % (allmotif[k][ind[i][k]],allstart[k][ind[i][k]],allstop[k][ind[i][k]],allEvalue[k][ind[i][k]]))
    f.write("\n")
    
f.close()  
