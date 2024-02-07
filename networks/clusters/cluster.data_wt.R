remove(list = ls()) 

NUM_CLUSTERS.REF <- 3

cluster_names <- c('CL1', 'CL2', 'CL3')
names(cluster_names) <- c('AML', 'Hybrid', 'Untreated')

library(sRACIPE)
source('./functions.racipe.wt.R')

# Load WT racipe models
#-----------------------
circuit_idx <- "0.07-32-0.75"  # best network 
#circuit_idx <- "0.08-60-0.3"  # largest network 
#circuit_idx <- "0.05-4-0.85"  # smallest network 

fname <- paste('../circuits.sim/circuit_simulated_', circuit_idx, '.rds', sep = '')
racipe <- readRDS(file = fname)

geneExpression <- assay(racipe,1) 
geneExpression <- log2(1+geneExpression)
means <- rowMeans(geneExpression)
sds <-  apply(geneExpression, 1, sd)

ndata <- normalize_by_wt_mean_and_sd(geneExpression=geneExpression, 
                                     means=means, 
                                     sds=sds) 
pca <- prcomp(t(ndata)) 

clusterCut <- obtain_clusterCut(ndata=ndata, num_clusters=NUM_CLUSTERS.REF)

sum(clusterCut==1)  
sum(clusterCut==2)  
sum(clusterCut==3)  

sum(clusterCut==1)*100/length(clusterCut)  
sum(clusterCut==2)*100/length(clusterCut)  
sum(clusterCut==3)*100/length(clusterCut)  


# Separate data into clusters 
#----------------------------
data_by_cluster <- obtain_data_by_cluster(data.sim=ndata, 
                                          clusterCut=clusterCut)
names(data_by_cluster)
