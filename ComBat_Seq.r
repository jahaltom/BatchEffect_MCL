library("sva") #Note this exercise requires sva (>= v3.36.0) which is only available for R (>= v4.x)
library("ggplot2")
library("gridExtra")
library("readr")
library(tidyr)
library(dplyr)
library(data.table)
library(edgeR)



#read metadata
metadata<-read_delim("AllCovidMetadata.txt","\t",escape_double=F,trim_ws=T, col_types = cols(.default = "c") )
#do this last
rownames(metadata)<-metadata$SampleID
dim(metadata)


#remove single batch rows
# metadata_reduced <- subset(metadata,duplicated(study_accession) | duplicated(study_accession, fromLast=TRUE))
# rownames(metadata_reduced)<-metadata_reduced$SampleID
# dim(metadata_reduced)



metadata=metadata[(metadata$'covid status-curated' == 'covid' | metadata$'covid status-curated' == 'non'),]


metadata_final <- metadata


counts=fread("SarsCov2_Studies_Counts.tsv",sep="\t", header=TRUE, stringsAsFactors=FALSE, showProgress=TRUE, nThread=30)
# limit counts to sample_name
#counts<-counts[ ,colnames(counts) %in% c('Gene_ID_ver',sample_names$V1), with=FALSE]
counts<-counts[ ,colnames(counts) %in% c('Gene_ID_ver',metadata_final$SampleID), with=FALSE]

#EB.chr12G3631, EB.chr4G15922, EB.chr7G54834
counts=counts[!grepl("EB.chr12G3631", counts$Gene_ID_ver),]
counts=counts[!grepl("EB.chr4G15922", counts$Gene_ID_ver),]
counts=counts[!grepl("EB.chr7G54834", counts$Gene_ID_ver),]


# filter low expressed genes
cpmdf<-cpm(counts %>% select(- c('Gene_ID_ver')))
# keep genes at least 10 cpm in at least 100 samples
tokeep<-rowSums(cpmdf >= 10) >= 100
#tokeep<-rowSums(cpmdf > 1) >= 10
# filter counts
counts<-counts[tokeep,]



#add gene names as rownames
rownames(counts)<-counts$Gene_ID_ver
gene_names<-counts$Gene_ID_ver
fwrite(list(gene_names),"genenames_order.txt")
counts$Gene_ID_ver <- NULL

#rearrange metadata rows to be identical as count
metadata_final <- metadata_final[match(names(counts), metadata_final$SampleID),]
rownames(metadata_final)<-metadata_final$SampleID
all(rownames(metadata_final) == colnames(counts))
fwrite(list(colnames(counts)),"samplenames_order.txt")



metadata_final$study_accession=metadata_final$study_accession %>% replace_na('MasaonAutopsy')
#running combat
#run combat seq
#adjusted <- ComBat_seq(counts, batch=metadata_final$BatchID, group=metadata_final$bio_group)
#adjusted <- ComBat_seq(counts, batch=metadata_final$BatchID, group=metadata_final$TissueType_details)
adjusted <- ComBat_seq(as.matrix(counts), batch=metadata_final$study_accession,group=metadata_final$'covid status-curated')


#AdjCounts to data table
adjusteddt <- as.data.table(adjusted)
adjusteddt$Gene_ID_ver<-gene_names
adjusteddt <- adjusteddt %>%  select(Gene_ID_ver, everything())
fwrite(adjusteddt, "SarsCov2_adjusted_counts.tsv", row.names=F, quote=FALSE, sep="\t")

#Reg Counts to data table
counts <- as.data.table(counts)
counts$Gene_ID_ver<-gene_names
counts <- counts %>%  select(Gene_ID_ver, everything())
fwrite(counts, "SarsCov2_Regular_counts.tsv", row.names=F, quote=FALSE, sep="\t")



save.image(file="adjdata.RData")
