library(clusterProfiler)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


df = read.csv('MCL', header=TRUE,sep='\t')


#Gather clusters >= 10
clustrsOrg=c()
for (i in 1:nrow(df)){
     #Loop through rows (clusters) df[row,col]   
    cluster= unlist(strsplit(df[i,2], split = ";"))
    
    #Remove EB and SarsCov2
    cluster=cluster[!cluster %in% grep(paste0("EB.chr", collapse = "|"), cluster, value = T)]
    cluster=cluster[!cluster %in% grep(paste0("SarsCov2", collapse = "|"), cluster, value = T)]
    
    if (length(cluster)>=10){

    clustrsOrg=c(clustrsOrg,list(cluster))}
}




#100 interations!
for (iteration in 1:100){    
    clusters=clustrsOrg
    #Randomize clusters with 100K gene swaps
    for (i in 1:100000){
        #Two random clusters
        c=sample(1:length(clusters), 2,replace=F)
        
        #Random genes in each cluster
        g1=sample(1:length(unlist(clusters[c[1]])), 1)
        g2=sample(1:length(unlist(clusters[c[2]])), 1)
        
        
        #Grab genes from each cluster
        rGeneC1=unlist(clusters[c[1]])[g1]  
        rGeneC2=unlist(clusters[c[2]])[g2]
           
            
        #Swap
        rc1=unlist(clusters[c[1]])
        rc1[g1]=rGeneC2
        
        rc2=unlist(clusters[c[2]])
        rc2[g2]=rGeneC1
        
        #Reinsert into origonal
        clusters[c[1]]=list(rc1)
        clusters[c[2]]=list(rc2)
    
    }
    
    
    
    
    
        
    #Go Enricment
    
    #For storing best p.adjust
    padj=c()
    
    for (i in 1:length(clusters)){
    
    

            gse=enrichGO(
            unlist(clusters[i]),
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
