

## Batch correction ComBatSeq

* ComBat_Seq.r: Takes in SarsCov2 counts (filteres to keep genes with at least 10 cpm in at least 100 samples) and assosisated metadata. Performs batch correction using study as "batch" and "group" as covid status. Outputs SarsCov2_adjusted_counts.tsv and SarsCov2_Regular_counts.tsv. 
```
metadata: All covid/non-covid data from AllCovid/CovidMetadata.xlsx.
counts: AllCovid/SarsCov2_Studies_Counts.tsv
```
 
## Gene-Gene correlation matrix and Markov Chain Clustering

## submitPearson.sh and submitSpearman.sh: Main workflow scripts. These do the folowing:

# gene-Gene correlation matrix

*  compute_pearson_sc.py and compute_spearman_sc.py: Take in SarsCov2_adjusted_counts.tsv or SarsCov2_Regular_counts.tsv from ComBat_Seq.r, and calculate the gene-gene pearson and spearman correlation matrix (respectfully). Outputs 3 column file (gene1IDVer   gene2IDVer  corr)

* For each matrix:
  * 3 files are generated from cutoffs (0.8, 0.85, 0.9) 
  * each resulting file is annotated by gene names (was gene ID Ver) using add_gene_names.py. Needs Gene_level_metadata.tsv. Outputs 3 column file (gene1Name   gene2Name  corr).

# Markov Chain Clustering

Each resulting file from above is ran through doMCL.py which generates the gene clusters in a 2 column file (ClusterID       Genes). 





## Go Enrichment with clusterProfiler

For each cluster, the SarsCov2 genes and EB genes are removed. Then clusteres with  >= 10 genes left are kept. The resulting gene clusteres are ran through "enrichGO" 
```
gse=enrichGO(
        cluster,
        OrgDb=organism,
        keyType = "SYMBOL",
        ont = "BP",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.2,
        minGSSize = 10,
        maxGSSize = 500,
        readable = FALSE,
        pool = FALSE
        )
```

For each cluster, sig results are exported. 

The mean of the best p.adjust for each cluster is recorded. 
