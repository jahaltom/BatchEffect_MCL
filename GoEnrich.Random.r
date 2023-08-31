library(clusterProfiler)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


df = read.csv('scorrs_thresh_0.9_renamed.tsv_mclout_inf_1.5', header=TRUE,sep='\t')


#Gather clusters of annoated genes
clustrsOrg=c()
#List of all annoated genes
annList=c()
for (i in 1:nrow(df)){
     #Loop through rows (clusters) df[row,col]   
    cluster= unlist(strsplit(df[i,2], split = ";"))
    
    if (length(cluster)>=10){
        #Remove EB and SarsCov2
        cluster=cluster[!cluster %in% grep(paste0("EB.chr", collapse = "|"), cluster, value = T)]
        cluster=cluster[!cluster %in% grep(paste0("SarsCov2", collapse = "|"), cluster, value = T)]

        clustrsOrg=c(clustrsOrg,list(cluster))
        annList=c(annList,cluster)}
}


#Remove empty clusters
clustrsOrg=clustrsOrg[lapply(clustrsOrg,length)>0]



#For 100 interations!
for (iteration in 1:100){  
    #Randomize list of all annoated genes
    annListR=sample(annList)
    clustrsOrgR=clustrsOrg  
    start=1
    for (c in 1:length(clustrsOrgR)){
            cLength=length(unlist(clustrsOrgR[c[1]]))
            clustrsOrgR[c[1]]=list(annListR[start:(cLength+start-1)])
            start=start+cLength
}   
        
    #Go Enricment
     
    #For storing best p.adjust
    padj=c()
    
    for (i in 1:length(clustrsOrgR)){
            gse=enrichGO(
            unlist(clustrsOrgR[i]),
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
            
            #Get min p-value for each cluster. If there are any. 
            if (length(gse$p.adjust) != 0){
                padj=c(padj,min(gse$p.adjust))}            
         }
        
        ########  
        #Average best p.adjust
        bestAvg=mean(padj)
        
        #Write output
        fileConn<-file("results/MCL_BestPValue.Random.txt","a")
        writeLines(as.character(bestAvg), fileConn)
        close(fileConn)

}
