#!/usr/bin/env Rscript

remove(list = ls()) 

NO_SAMPLES <- 2 #1000
SAMPLE_SIZE <- 10000

NO.TOP.CIRCUITS <- 3 #10
PERCENT_REDUCTION <- 90

library(sRACIPE)

source('../../networks/functions.R')

outdir <- './results/'
dir.create(outdir)

circuit_metrics.sim <- read.csv(file = '../../networks/results/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)
circuit_metrics.sim <- circuit_metrics.sim[1:NO.TOP.CIRCUITS, ]

flexibility.df <- as.data.frame(matrix(nrow = NO.TOP.CIRCUITS, ncol=NO_SAMPLES))
dim(flexibility.df)
rownames(flexibility.df) <- rownames(circuit_metrics.sim)

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

rv = foreach(circuit_idx = rownames(circuit_metrics.sim), .combine = 'rbind', .inorder = TRUE) %dopar% { 
  library(sRACIPE)
  racipe.all <- readRDS(file = paste('../networks/circuits.sim/circuit_simulated_', circuit_idx, '.rds', sep = ''))
  flex.v <- sapply(seq(1:NO_SAMPLES), function(x){
      racipe <- racipe.all[, sample(ncol(racipe.all), SAMPLE_SIZE)] 
      racipe <- sRACIPE::sracipeNormalize(racipe) 
      # KD by subsetting:
      racipe.kd <- sRACIPE::sracipeKnockDown(racipe, plotToFile = FALSE,
                                           plotBarPlot = FALSE, #TRUE, 
                                           plotHeatmap = FALSE, 
                                           reduceProduction = (100-PERCENT_REDUCTION)
    )
    avg.dist <- calDistance(racipe.kd)
    return(avg.dist)
  })
  write.csv(flex.v, file = paste(outdir, "./flexibility-", circuit_idx,".csv", sep = ''),
            row.names = T, quote = F)
  flex.v
}

rownames(rv) <- rownames(circuit_metrics.sim) 
colnames(rv) <- paste('V', seq(1, NO_SAMPLES), sep = '')
write.csv(rv, file = paste0(outdir, "./flexibility.all.csv"), row.names = T, quote = F) 

# Release resources 
stopCluster(cl)

