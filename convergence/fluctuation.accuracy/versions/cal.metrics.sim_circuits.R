#!/usr/bin/env Rscript

remove(list = ls()) 

# PARAMETERS
# NUM_CLUSTERS: user supplied information about 
# how many clusters in the data
NUM_CLUSTERS <- 2

library(igraph)
library(NetAct)
library(sRACIPE)
source('./functions.R')
source('./heatmapSimilarity.updated.R')
library(dplyr)

outdir.hS <- './circuits.hS/'
dir.create(outdir.hS) 

outdir <- './results/'
dir.create(outdir)

# Load circuit summary
#-------------------- 
circuit_metrics <- read.csv(file = './circuits/summary.circuits.csv', row.names = 1)

# select section of the metrics specific to the simulated circuits:
circuit_metrics.sim <- circuit_metrics[!circuit_metrics$DupStatus,]

# Load circuits
circuits <- readRDS('./circuits/circuits.uniq.rds')
names(circuits)
length(names(circuits))
fname.list <- list() 
for(idx in names(circuits)){
   print(idx)
   fname.list[[idx]] <- paste('./circuits.sim/circuit_simulated_', idx, '.rds', sep = '')
}
names(fname.list) <- names(circuits)


# select TOP n TFs from each of the THREE methods: NetAct, MARINa, and RI
#------------------------------------------------------------------------
coreTFs.list <- readRDS('../tfSets/data/coreTFs.rds')
names(coreTFs.list)

# Construct possible circuit set for the TFs based on each TF-target DB in the TF-target DB set
#---------------------------------------------------------------------------------------------- 
targetDB.list <- readRDS('../databases/targetDB.list.rds')
names(targetDB.list)


# Load ExpressionSet object and DE results ---- START 
# Load brain array expression data as an ExpressionSet object
#------------------------------------------------------------
fname.eset.brain_array <- '../data.tfs/eset.brain_array.rda'
load(fname.eset.brain_array) # loads object eset.brain_array  
# Load DE results 
#---------------- 
fname.de.results <- '../data.tfs/de.results.rda'
load(fname.de.results) # loads object de.results 
# Load ExpressionSet object and DE results ---- END 



# Load cluster cut for reference data:
#----------------------------------------
#clusterCut.REF <- clusterCut
#write.csv(clusterCut.REF, file = paste0('./circuits/', 'clusterCut.REF.csv'), quote = F, row.names = T) 

tmp <- read.csv(file = paste0('./circuits/', 'clusterCut.REF.csv'), row.names = 1)
clusterCut.REF <- as.integer(tmp$x)
names(clusterCut.REF) <- as.character(rownames(tmp))


# Calculate heatmap similarity between data.REF and data.sim 
#------------------------------------------------------------
for(circuit_idx in names(circuits)){
   #circuit_idx <- names(circuits)[1]
   #print(circuit_idx)
   #break()
   # obtain TF activities using NetAct:
   fr <- strsplit(circuit_idx, split = '-')[[1]][1]
   top.TFs.count <- strsplit(circuit_idx, split = '-')[[1]][2]  
   
   #coreTFs <- coreTFs.list[[top.TFs.count]] 
   coreTFs <- coreTFs.list[[fr]][[top.TFs.count]][['COMB']]
   targetDB = targetDB.list[[fr]] 
   length(coreTFs)
   coreTFs <- intersect(coreTFs, names(targetDB))
   length(coreTFs) 
   a = TF_Activity(tfs = coreTFs,
                   GSDB = targetDB,  
                   eset = eset.brain_array,
                   DErslt = de.results  #DErslt=de.results$Overall
   )
   data.REF <- a$all_activities
   
   # obtain simulated data:
   racipe <- readRDS(file = paste('./circuits.sim/circuit_simulated_', 
                                        circuit_idx, '.rds', sep = '')) 
   racipe <- sracipeNormalize(racipe)  
   data.sim <- assay(racipe,1)  
   
   #break()
   # calculate similarity between activitites and racipe simulation data:
   hS <- sracipeHeatmapSimilarity(dataReference = data.REF,
                                  dataSimulation = data.sim, 
                                  returnData = T, 
                                  #nClusters = NUM_CLUSTERS, 
                                  clusterCut = clusterCut.REF)  
   
   # add similarity values to the metrics data frame:
   circuit_metrics.sim[circuit_idx, "ClusterA"] <- hS$simulated.cluster.freq[2] 
   circuit_metrics.sim[circuit_idx, "ClusterB"] <- hS$simulated.cluster.freq[3] 
   circuit_metrics.sim[circuit_idx, "Accuracy"] <- (hS$simulated.cluster.freq[2] + hS$simulated.cluster.freq[3])
   circuit_metrics.sim[circuit_idx, "KLdist"] <- hS$KL 
   circuit_metrics.sim[circuit_idx, "ClusterA2"] <- hS$cluster.similarity[2] 
   circuit_metrics.sim[circuit_idx, "ClusterB2"] <- hS$cluster.similarity[3] 
   circuit_metrics.sim[circuit_idx, "AvgDist"] <- hS$AvgDist 
   
   # save hS object to an output file:
   fname.hS <- paste(outdir.hS, 'hS_', circuit_idx, '.rds', sep = '')
   saveRDS(hS, file = fname.hS)
   #break()
}

write.csv(circuit_metrics.sim, file = paste0(outdir, "./summary.circuits.sim.csv"),
          row.names = T, quote = F) 

