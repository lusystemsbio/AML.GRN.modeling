#!/usr/bin/env Rscript

remove(list = ls()) 

NO_SAMPLES <- 2 #1000
SAMPLE_SIZE <- 10000

NO.TOP.CIRCUITS <- 10
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


for(circuit_idx in rownames(circuit_metrics.sim)){
  print(circuit_idx) 
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
  flexibility.df[circuit_idx, ] <- flex.v
  #break()
}

write.csv(flexibility.df, file = paste0(outdir, "./flexibility.seq.csv"), row.names = T, quote = F)




