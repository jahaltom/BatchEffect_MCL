import sys
import pandas as pd
import numpy as np




df=pd.read_csv(sys.argv[1],sep='\t',skiprows=0)
gene_idVer=df['Gene_ID_ver'].tolist()
df = df.set_index('Gene_ID_ver')
npdf=df.to_numpy()
corrs=np.corrcoef(npdf)


headr='gene1\tgene2\tspearman_corr'
print(headr)
for idx in range(len(gene_idVer)):
    currgene=gene_idVer[idx]
    for j in range(idx,len(gene_idVer)):
        currcorr=round(corrs[idx][j],4)      
        refgene=gene_idVer[j]
        print('\t'.join([currgene,refgene,str(currcorr)]))






