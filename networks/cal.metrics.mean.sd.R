#!/usr/bin/env Rscript
remove(list = ls()) 

# PARAMETERS

# NUM_CLUSTERS: user supplied information about 
# how many clusters in the data
NUM_CLUSTERS <- 2
SAMPLE_SIZE <- 2000

library(sRACIPE)
source('./functions.R')
source('./heatmapSimilarity.updated.R')

outdir <- './circuits/'
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

# Load cluster cut for reference data:
#----------------------------------------
tmp <- read.csv(file = paste0('./circuits/', 'clusterCut.REF.csv'), row.names = 1)
clusterCut.REF <- as.integer(tmp$x)
names(clusterCut.REF) <- as.character(rownames(tmp))


# Load TF activities 
#--------------------
fname.TFactivities <- './circuits/circuit_TFactivities.rds'
TF_activities <- readRDS(fname.TFactivities) 

# Calculate heatmap similarity between data.REF and data.sim 
#------------------------------------------------------------
for(circuit_idx in names(circuits)){
   data.REF <-  TF_activities[[circuit_idx]]
   # obtain simulated data:
   racipe <- readRDS(file = paste('./circuits.sim/circuit_simulated_', 
                                        circuit_idx, '.rds', sep = '')) 
   racipe <- sracipeNormalize(racipe)  
   data.sim <- assay(racipe,1)  

   AvgDist <- cal.avgDist.bySample(data.REF, data.sim, clusterCut.REF, SAMPLE_SIZE) 
   circuit_metrics.sim[circuit_idx, "AvgDist.mean"] <- mean(AvgDist)
   circuit_metrics.sim[circuit_idx, "AvgDist.sd"] <- sd(AvgDist)

   #break()
}

write.csv(circuit_metrics.sim, file = paste0(outdir, './summary.circuits.sim.mean_sd.csv'),
          row.names = T, quote = F) 

