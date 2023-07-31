import sys

genemd="Gene_level_metadata.tsv"



genedict={}
with open(genemd) as fi:
    for l in fi:
        tmp=l.split('\t')
        genedict[tmp[22]]=tmp[1]

ind=0
with open("scorrs_thresh_0.8.tsv") as fi:
    for line in fi:
        l=line.strip()
        if ind==0:
            print(l)
            ind+=1
            continue
        tmp=l.split('\t')
        tmp[0]=genedict[tmp[0]]
        tmp[1]=genedict[tmp[1]]
        print('\t'.join(tmp))
