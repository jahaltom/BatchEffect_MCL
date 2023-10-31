import pandas as pd






#Read in cluster file with all annoated and EB genes.
clusterFile=pd.read_csv('scorrs_thresh_0.85_renamed.tsv_mclout_inf_1.6',sep='\t')
clusters= clusterFile["ClusterID"].tolist()     


                
#To store GoTerms for all clusters
go=[]
#Loop through clusters
for i in clusters:
    try:
        #grab top 10 GoTerms. 
        goFile=pd.read_csv('results/scorrs_thresh_0.85_renamed.tsv_mclout_inf_1.6_Cluster_'+str(i) +'_GoEnrichmentResults',sep='\t').head(10)
        #Gather cluster info (ClusterID, Genes) and add to GoTerm info for the same cluster. 
        cluster=clusterFile[clusterFile["ClusterID"]==i]         
        goFile["Cluster"]=cluster["ClusterID"].tolist()[0]
        goFile["Genes"]=cluster["Genes"].tolist()[0]
        #Append to master list
        go.append(goFile)
    except:
        #Gather cluster info (ClusterID, Genes) and add to GoTerm info for the same cluster. 
        cluster=clusterFile[clusterFile["ClusterID"]==i]         
        cluster = cluster.rename(columns={'ClusterID': 'Cluster'})
        #Append to master list
        go.append(cluster)
        
    
    
go=pd.concat(go)

go.to_csv("out.tsv",sep='\t')
