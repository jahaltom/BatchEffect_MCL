import sys
import pandas as pd
import numpy as np
from scipy.stats import spearmanr



#print('reading data')
df=pd.read_csv("SarsCov2_adjusted_counts.tsv",sep='\t',skiprows=0)
gene_names=df['Gene_ID_ver'].tolist()
#print('GN',gene_names)
df = df.set_index('Gene_ID_ver')
#print(df)
npdf=df.to_numpy()

res=spearmanr(npdf,axis=1)
#print(res)
corrs=res.correlation
pvals=res.pvalue

headr='gene1\tgene2\tspearman_corr\tpval'
print(headr)
for idx in range(len(gene_names)):
    currgene=gene_names[idx]
    #currcorr=corrs[idx][idx:]
    for j in range(idx,len(gene_names)):
        currcorr=round(corrs[idx][j],4)
        currpval=round(pvals[idx][j],4)
        refgene=gene_names[j]
        print('\t'.join([currgene,refgene,str(currcorr), str(currpval)]))
    #print(currgene,currcorr)
#write to file
#np.savetxt('npput.tsv',res.correlation,delimiter='\t')
#print('done')



