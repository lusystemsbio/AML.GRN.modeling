#!/usr/bin/env Rscript
# Resources: 
# 1. http://r-statistics.co/Linear-Regression.html 

remove(list = ls())

circuit_metrics.sim <- read.csv(file = './results/summary.circuits.sim.sigVal.csv', row.names = 1)

source('../../networks/functions.R')

outdir <- './results/'

# Three way sorting - accuracy and flexibility distance
#======================================================
circuit_metrics.sim.sorted <- sortByTwoIndices.acc_and_flex(circuit_metrics.sim) 
write.csv(circuit_metrics.sim.sorted, 
          file = './results/summary.circuits.sim.sortedByAcc_flex.csv',  
          row.names = T, quote = F)  