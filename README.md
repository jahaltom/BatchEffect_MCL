
### Predicting the function of human evidence based genes involved in SarsCov2.

Human transcriptome(GencodeV36), SARS-COV-2 transcriptome (ASM985889v3), and human evidence based (EB) gene transcripts were used for transcriptome mapping (https://github.com/jahaltom/COVID-19-Quantification/tree/main) using 66 human SarsCov2 RNA-Seq studies (3,413 samples). To increase mapping accuracy, the human genome along with viral decoys and spike-ins from the Genomic Data Commons (GRCh38.d1.vd1) were used as decoys for Salmon. 

This massive amount of RNA-Seq expression data can be mined to look for EB co-expression with annotated genes and this can shed insights into EB gene function. The co-expression can be elucidated by performing gene-gene correlation analysis, in this case we performed spearman and pearson. We did this with both raw counts and ComBatSeq batch corrected counts. The subsequent correlation matrices can then be put through Markov Chain Clustering (MCL) to find clusters of genes with correlated expression across the whole transcriptome. The annotated genes, tied to the EB genes, in the MCL clusters can be ran through GoEnrichment analysis to determine possible EB gene functions. 

This method is further established by taking the exact same clusters and randomly shuffling the genes throughout the clusters ( keeping the number of clusters and number of genes within each cluster the same) for 100 iterations. Performing GoEnrichment analysis on these randomized clusters reviled significantly less enrichment. 



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

# Experimental GoEnrich.r
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

# Random GoEnrich.Random.r

Same as above except there is no sig results exported, and the clusteres are randomized after the ">=10"  step (prior to "enrichGO"). 

* Randomization: Two random cluteres are pulled out and a single random gene from each cluster are swapped. This is done 100K times. 

The mean of the best p.adjust for each cluster is recorded.

This is done for 100 iterations. The 100 means will be ploted along with the Experimental.






























