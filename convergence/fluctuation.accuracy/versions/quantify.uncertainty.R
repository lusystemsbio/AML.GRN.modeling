#!/usr/bin/env Rscript 
rm(list = ls()) 

NO.TOP.CIRCUITS <- 10

datadir <- '../rankings.allCircuits/results/'
circuit_metrics.sim <- read.csv(file = paste0(datadir , "./summary.circuits.sim.sortedByAcc_flex.csv"), row.names = 1)


circuit_metrics.sim <- circuit_metrics.sim[1:NO.TOP.CIRCUITS, ]






