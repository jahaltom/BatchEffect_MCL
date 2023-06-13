import pandas as pd





#Read in cluster #s that contain EB(s)
clusters=pd.read_csv('EB_Clusters',sep='\t',header=None)  
clusters=clusters[0].tolist()

#Read in cluster file with all annoated and EB genes.
clusterFile=pd.read_csv('../scorrs_thresh_0.8_renamed.tsv_mclout_inf_1.9',sep='\t')
                        
#To store GoTerms for all clusters
go=[]
#Loop through clusters
for i in clusters:
    try:
        #grab top 10 GoTerms. 
        goFile=pd.read_csv('scorrs_thresh_0.8_renamed.tsv_mclout_inf_1.9_Cluster_'+str(i) +'_GoEnrichmentResults',sep='\t').head(10)
        #Gather cluster info (ClusterID, Genes) and add to GoTerm info for the same cluster. 
        cluster=clusterFile[clusterFile["ClusterID"]==i]         
        goFile["Cluster"]=cluster["ClusterID"].tolist()[0]
        goFile["Genes"]=cluster["Genes"].tolist()[0]
        #Append to master list
        go.append(goFile)
    except:
        pass
    
    
go=pd.concat(go)

go.to_csv("Pearson0.8_ComBatSeq_GoEnrichment.tsv",mode="w", header=True,index=False,sep="\t")
