
# Predicting the function of human evidence based genes involved in SarsCov2.

Human transcriptome(GencodeV36), SARS-COV-2 transcriptome (ASM985889v3), and human evidence based (EB) gene transcripts were used for transcriptome mapping (https://github.com/jahaltom/COVID-19-Quantification/tree/main) using 66 human SarsCov2 RNA-Seq studies (3,413 samples). To increase mapping accuracy, the human genome along with viral decoys and spike-ins from the Genomic Data Commons (GRCh38.d1.vd1) were used as decoys for Salmon. 

This massive amount of RNA-Seq expression data can be mined to look for EB co-expression with annotated genes and this can shed insights into EB gene function. The co-expression can be elucidated by performing gene-gene correlation analysis, in this case we performed spearman and pearson. We did this with both raw counts and ComBatSeq batch corrected counts. The subsequent correlation matrices can then be put through Markov Chain Clustering (MCL) to find clusters of genes with correlated expression across the whole transcriptome. The annotated genes, tied to the EB genes, in the MCL clusters can be ran through GoEnrichment analysis to determine possible EB gene functions. 

This method is further established by taking the exact same clusters and randomly shuffling the genes throughout the clusters ( keeping the number of clusters and number of genes within each cluster the same) for 100 iterations. Performing GoEnrichment analysis on these randomized clusters reviled significantly less enrichment. 



## Batch correction ComBatSeq

* ComBat_Seq.r: Takes in SarsCov2 counts "SarsCov2_Studies_Counts.tsv" and filteres to keep genes with at least 10 cpm in at least 100 samples). Also takes assosisated metadata AllCovidMetadata.txt. Performs batch correction using study as "batch" and "group" as covid status. Outputs SarsCov2_adjusted_counts.tsv and SarsCov2_Regular_counts.tsv. 

ComBat_Seq.r input:
```
metadata: All covid/non-covid data from AllCovid/CovidMetadata.xlsx.
counts: AllCovid/SarsCov2_Studies_Counts.tsv
```
 
## Gene-Gene correlation matrix and Markov Chain Clustering
## submitPearson.sh and submitSpearman.sh: Main workflow scripts. These do the folowing:

each combonation (Person/spearman regular/adjusted counts = 4 total)  need own dir. 
### gene-Gene correlation matrix

*  compute_pearson_sc.py and compute_spearman_sc.py: Take in SarsCov2_adjusted_counts.tsv or SarsCov2_Regular_counts.tsv from ComBat_Seq.r, and calculate the gene-gene pearson and spearman correlation matrix (respectfully). Outputs 3 column  (gene1IDVer   gene2IDVer  corr) file:  scorrs.tsv.

* For each matrix:
  * 3 files are generated from corr cutoffs (0.8, 0.85, 0.9): scorrs_thresh_0.9.tsv, scorrs_thresh_0.8.tsv, scorrs_thresh_0.85.tsv
  * each resulting file is annotated by gene names (was gene ID Ver) using add_gene_names.py. Needs Gene_level_metadata.tsv. Outputs 3 column (gene1Name   gene2Name  corr) files: scorrs_thresh_0.9_renamed.tsv, scorrs_thresh_0.8_renamed.tsv, scorrs_thresh_0.85_renamed.tsv.

### Markov Chain Clustering

Each resulting file from above is ran through doMCL.py which generates the gene clusters in a 2 column (ClusterID       Genes) files:  scorrs_thresh_0.9_renamed.tsv_mclout_inf_1.5, scorrs_thresh_0.8_renamed.tsv_mclout_inf_1.5, scorrs_thresh_0.85_renamed.tsv_mclout_inf_1.5.





## Go Enrichment with clusterProfiler

### Experimental: GoEnrich.r
For each cluster, the SarsCov2 genes and EB genes are removed. Then the clusters are filtered to contain at least 10 known genes. The resulting gene clusteres are ran through "enrichGO" 
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

### Random: GoEnrich.Random.r

Same as above except there is no sig results exported, and the clusteres are randomized. 

* Randomization: All genes from all clusters are extracted and randomized into a list. Then each cluster is repopulated with the randomized gene list, keeping the cluster sizes the same as the original.

The mean of the best p.adjust for each cluster is recorded.

This is done for 100 iterations. The 100 means will be ploted along with the Experimental.


### Run all with submit.sh and submit.Random.sh


```
mkdir results


cat list | while read i;do
    cat GoEnrich.r | sed "s/MCL/$i/g" > $i.GoEnrich.r
    cat submit.sh | sed "s/GoEnrich.r/$i.GoEnrich.r/g" > $i.submit.sh
    sbatch $i.submit.sh
done


cat list | while read i;do
    cat GoEnrich.Random.r | sed "s/MCL/$i/g" > $i.GoEnrich.Random.r
    cat submit.Random.sh | sed "s/GoEnrich.Random.r/$i.GoEnrich.Random.r/g" > $i.submit.Random.sh
    sbatch $i.submit.Random.sh
done
```











## R PLots
```
#For each output.Random.txt and output.txt
ls *txt* | cat | sed 's/_output.Random.txt//g' | sed 's/_output.txt//g' | sort  | uniq > files





cat files | while read i;do
    cat "$i"_output.Random.txt | grep Mean | awk '{print $NF}' > rand
    cat "$i"_output.txt | grep Mean | awk '{print $NF}' > exp
    sed "s/FILE/$i/g" Plots.r > $i.Plots.r
    Rscript $i.Plots.r
done
```



## GoTerm analysis
```
#Gather all clusters that contain EB genes
cat scorrs_thresh_0.8_renamed.tsv_mclout_inf_1.9 | grep "EB.chr" | awk '{print $1}' > EB_Clusters

Rscript ClusterGo.r

```




