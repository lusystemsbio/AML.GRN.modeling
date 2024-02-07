#!/usr/bin/env Rscript
remove(list = ls()) 

NO_SAMPLES <- 10
SAMPLE_SIZE <- 10000
NO_TOP_CIRCUITS <- 10
PERCENT_REDUCTION <- 95

library(sRACIPE)
source('../../networks/functions.R')
source('../../networks/heatmapSimilarity.updated.R') 
source('./functions.R')

outdir <- './results/'
dir.create(outdir)

# Load top circuits
#------------------
#circuits <- readRDS('../../../networks/circuits/circuits.uniq.rds')
#circuit_metrics.sim.sorted <- read.csv('../../../networks/results/summary.circuits.sim.sortedByAcc_flex.csv', row.names=1)
circuit_metrics.sim.sorted <- read.csv('../rankings.mean_acc/results/summary.circuits.sim.sortedByAcc_flex.csv', row.names=1)
#circuits <- rownames(circuit_metrics.sim.sorted)[2:(NO_TOP_CIRCUITS+1)]
circuits <- rownames(circuit_metrics.sim.sorted)[1:(NO_TOP_CIRCUITS)]
names(circuits) <- circuits

# Load circuit summary
#-------------------- 
circuit_metrics <- read.csv(file = '../../networks/circuits/summary.circuits.csv', row.names = 1)
# select section of the metrics specific to the simulated circuits:
circuit_metrics.sim <- circuit_metrics[circuits, ]

#names(circuits)
# length(names(circuits))
# fname.list <- list() 
# for(idx in names(circuits)){
#    #print(idx)
#    fname.list[[idx]] <- paste('../../circuits.sim/circuit_simulated_', idx, '.rds', sep = '')
# }
# names(fname.list) <- names(circuits)

# Load cluster cut for reference data:
#----------------------------------------
tmp <- read.csv(file = paste0('../../networks/circuits/', 'clusterCut.REF.csv'), row.names = 1)
clusterCut.REF <- as.integer(tmp$x)
names(clusterCut.REF) <- as.character(rownames(tmp))

# Load TF activities 
#--------------------
fname.TFactivities <- '../../networks/circuits/circuit_TFactivities.rds'
TF_activities <- readRDS(fname.TFactivities) 

# Calculate accuracy and average distance for 10 samples 
# For each sample, calculate these metrics for all circuits
#---------------------------------------------------------- 
circuit_metrics.sim.list <- list()
for(SAMPLE_NO in 1:NO_SAMPLES){
   sample_name <- paste0('S', SAMPLE_NO) 
   start.idx <- SAMPLE_SIZE * (SAMPLE_NO-1) + 1
   end.idx <- SAMPLE_SIZE * SAMPLE_NO
   # print(start.idx)
   # print(end.idx) 
   
   # Create data structure to save metrics for the current samples:
   #--------------------------------------------------------------
   # accuracy.avgDist <- as.data.frame(matrix(nrow = dim(circuit_metrics.sim)[1], ncol = 0) )
   # accuracy.avgDist$SimIdx <- circuit_metrics.sim$SimIdx
   # rownames(accuracy.avgDist) <- rownames(circuit_metrics.sim) 
   
   accuracy.avgDist <- circuit_metrics.sim 
   
   # Calculate metricx (Accuracy and AvgDist) for the current sample:
   #----------------------------------------------------------------
   circuit_idx <- names(circuits)[1]
   for(circuit_idx in names(circuits)){
      print(circuit_idx)
      data.REF <-  TF_activities[[circuit_idx]]
      # obtain simulated data:
      racipe <- readRDS(file = paste('../../networks/circuits.sim/circuit_simulated_', 
                                     circuit_idx, '.rds', sep = ''))  
      racipe <- sracipeNormalize(racipe)  
      racipe <- racipe[, start.idx:end.idx]  
      data.sim <- assay(racipe,1)  
      #data.sim.sub <- data.sim[, start.idx:end.idx] 
      #break() 
      
      # calculate similarity between activitites and racipe simulation data:
      hS <- sracipeHeatmapSimilarity(dataReference = data.REF,
                                     dataSimulation = data.sim, 
                                     returnData = T, 
                                     #nClusters = NUM_CLUSTERS, 
                                     clusterCut = clusterCut.REF)
      
      accuracy.avgDist[circuit_idx, "ClusterA"] <- hS$simulated.cluster.freq[2] 
      accuracy.avgDist[circuit_idx, "ClusterB"] <- hS$simulated.cluster.freq[3]
      accuracy.avgDist[circuit_idx,"KLdist"] <- hS$KL 
      accuracy.avgDist[circuit_idx,"ClusterA2"] <- hS$cluster.similarity[2] 
      accuracy.avgDist[circuit_idx,"ClusterB2"] <- hS$cluster.similarity[3]
      accuracy.avgDist[circuit_idx, "Accuracy"] <- (hS$simulated.cluster.freq[2] + hS$simulated.cluster.freq[3])
      accuracy.avgDist[circuit_idx, "AvgDist"] <- hS$AvgDist  
      
      # calculate circuit flexibility by KD subsetting  
      racipe.kd <- sRACIPE::sracipeKnockDown(racipe, 
                                             plotToFile = FALSE,
                                             plotBarPlot = FALSE, #TRUE,
                                             plotHeatmap = FALSE,
                                             reduceProduction = (100-PERCENT_REDUCTION)
      )
      flexibility_as_avg.dist <- calDistance(racipe.kd)
      accuracy.avgDist[circuit_idx, "flexibility"] <- flexibility_as_avg.dist

      #break()
   }
   circuit_metrics.sim.list[[sample_name]] <- accuracy.avgDist
   #break()
}

names(circuit_metrics.sim.list)
saveRDS(circuit_metrics.sim.list, file = paste0(outdir,  './circuit_metrics.sim.list.rds'))
