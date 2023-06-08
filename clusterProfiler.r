library(clusterProfiler)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


df = read.csv('scorrs_thresh_0.8_renamed.tsv_mclout_inf_1.9', header=TRUE,sep='\t')


#For storing best p.adjust
padj=c()
for (i in 1:nrow(df)){
     #Loop through rows (clusters) df[row,col]   
    cluster= unlist(strsplit(df[i,2], split = ";"))



    gse=enrichGO(
    cluster,
    OrgDb=organism,
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE
    )
    
    #Get min p-value for each cluster
    
    padj=c(padj,min(gse$p.adjust))
    
 }



#Average best p.adjust
best=mean(padj)

print(best,nrow(df))
