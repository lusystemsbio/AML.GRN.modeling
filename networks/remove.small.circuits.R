#!/usr/bin/env Rscript
# Resources: 
# 1. http://r-statistics.co/Linear-Regression.html 

remove(list = ls())
MINIMUM.CIRCUIT.SIZE <- 15
outdir <- './results/'

circuit_metrics.sim <- read.csv(file = './results/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)

circuit_metrics.sim.tmp <- circuit_metrics.sim[circuit_metrics.sim$Nodes>=MINIMUM.CIRCUIT.SIZE, ]

write.csv(format(circuit_metrics.sim.tmp, digits = 2), 
          file = './results/summary.circuits.sim.sortedByAcc_flex.nodes-15orMore.csv',  
          row.names = T, quote = F) 
