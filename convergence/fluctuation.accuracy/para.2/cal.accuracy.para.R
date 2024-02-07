#!/usr/bin/env Rscript

remove(list = ls()) 

# PARAMETERS
# NUM_CLUSTERS: user supplied information about 
# how many clusters in the data
NUM_CLUSTERS <- 2
NO.TOP.CIRCUITS <- 10
NO.REPEATS <- 1000
hS.PERMUTATIONS <- 1000

library(igraph)
library(NetAct)
library(sRACIPE)
source('../../../networks/functions.R')
source('../../../networks/heatmapSimilarity.updated.R')
library(dplyr)

outdir <- './results/'
dir.create(outdir)

# Load circuit summary
#-------------------- 
datadir <- '../../rankingsBypercentile/results/'
circuit_metrics.sim <- read.csv(file = paste0(datadir , "./summary.circuits.sim.sortedByAcc_flex.csv"), row.names = 1)
circuit_metrics.sim <- circuit_metrics.sim[1:NO.TOP.CIRCUITS, ]


# select TOP n TFs from each of the THREE methods: NetAct, MARINa, and RI
#------------------------------------------------------------------------
coreTFs.list <- readRDS('../../../tfSets/data/coreTFs.rds')
names(coreTFs.list)

# Construct possible circuit set for the TFs based on each TF-target DB in the TF-target DB set
#---------------------------------------------------------------------------------------------- 
targetDB.list <- readRDS('../../../databases/targetDB.list.rds')
names(targetDB.list)


# Load ExpressionSet object and DE results ---- START 
# Load brain array expression data as an ExpressionSet object
#------------------------------------------------------------
fname.eset.brain_array <- '../../../data.tfs/eset.brain_array.rda'
load(fname.eset.brain_array) # loads object eset.brain_array  
# Load DE results 
#---------------- 
fname.de.results <- '../../../data.tfs/de.results.rda'
load(fname.de.results) # loads object de.results 
# Load ExpressionSet object and DE results ---- END 


# Load cluster cut for reference data:
#----------------------------------------
#clusterCut.REF <- clusterCut
#write.csv(clusterCut.REF, file = paste0('./circuits/', 'clusterCut.REF.csv'), quote = F, row.names = T) 

tmp <- read.csv(file = paste0('../../../networks/circuits/', 'clusterCut.REF.csv'), row.names = 1)
clusterCut.REF <- as.integer(tmp$x)
names(clusterCut.REF) <- as.character(rownames(tmp))

# accuracy.df <- as.data.frame(matrix(nrow = nrow(circuit_metrics.sim), ncol = NO.REPEATS))
# rownames(accuracy.df) <- rownames(circuit_metrics.sim)

# Allocate resources for parallel computing
#-------------------------------------------#
library(doParallel)
library(foreach)
NO_AVAILABLE_CORES <- detectCores() 
print(NO_AVAILABLE_CORES)
cl <- makeCluster(NO_AVAILABLE_CORES)
#cl <- makeCluster(NO_AVAILABLE_CORES, setup_strategy = "sequential")
class(cl)
registerDoParallel(cl)
getDoParWorkers()

set.seed(1)
# Calculate heatmap similarity between data.REF and data.sim 
#------------------------------------------------------------
# for(circuit_idx in rownames(circuit_metrics.sim)){
rv = foreach(circuit_idx = rownames(circuit_metrics.sim), .combine = 'rbind', .inorder = TRUE) %dopar% { 
   library(sRACIPE)
   library(NetAct)
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
   racipe <- readRDS(file = paste('../../../networks/circuits.sim/circuit_simulated_', 
                                        circuit_idx, '.rds', sep = '')) 
   racipe <- sracipeNormalize(racipe)  
   data.sim <- assay(racipe,1)  
   
   # calculate similarity between activitites and racipe simulation data: 
   acc.v <- sapply(seq(1,NO.REPEATS),   
            function(x){ 
                  hS <- sracipeHeatmapSimilarity(dataReference = data.REF,
                                                  dataSimulation = data.sim, 
                                                  returnData = T, 
                                                  #nClusters = NUM_CLUSTERS, 
                                                  clusterCut = clusterCut.REF, 
                                                  permutations = hS.PERMUTATIONS)
                  return(hS$simulated.cluster.freq[2] + hS$simulated.cluster.freq[3])
          }
   )
   write.csv(acc.v, file = paste(outdir, "./accuracy-", circuit_idx,".csv", sep = ''),
             row.names = T, quote = F)
   acc.v
}
class(rv)
rownames(rv) <- rownames(circuit_metrics.sim) 
colnames(rv) <- paste('V', seq(1,NO.REPEATS), sep = '')
write.csv(rv, file = paste0(outdir, "./accuracy.all.csv"), row.names = T, quote = F) 

# Release resources for parallel computation 
#---------------------------------------------#
stopCluster(cl)
