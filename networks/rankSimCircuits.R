#!/usr/bin/env Rscript
# Resources: 
# 1. http://r-statistics.co/Linear-Regression.html 

remove(list = ls())

circuit_metrics.sim <- read.csv(file = './results/summary.circuits.sim.sigVal.csv', row.names = 1)

source('./functions.R')

outdir <- './results/'

# Four way sorting: 
#=================
#  (1) first by Accuracy and index the rows 
#  (2) second by AvgDist and index the rows 
#  (3) thrid by flexibility and index the rows
#  (4) add the three indices and sort by the summed indices for the fourth time
#circuit_metrics.sim <- sortByTwoIndices(circuit_metrics.sim) 
circuit_metrics.sim.sorted <- sortByThreeIndices(circuit_metrics.sim) 

sum(rownames(circuit_metrics.sim.sorted) =="0.09-50-0.75")


# Save the sorted circuits
#-------------------------
write.csv(circuit_metrics.sim.sorted, 
          file = './results/summary.circuits.sim.sorted.csv',  
          row.names = T, quote = F) 

# Three way sorting - accuracy and avg distance
#==============================================
circuit_metrics.sim.sorted <- sortByTwoIndices(circuit_metrics.sim) 
write.csv(circuit_metrics.sim.sorted, 
          file = './results/summary.circuits.sim.sortedByAcc_avgDist.csv',  
          row.names = T, quote = F) 

# Three way sorting - accuracy and flexibility distance
#======================================================
circuit_metrics.sim.sorted <- sortByTwoIndices.acc_and_flex(circuit_metrics.sim) 
write.csv(circuit_metrics.sim.sorted, 
          file = './results/summary.circuits.sim.sortedByAcc_flex.csv',  
          row.names = T, quote = F) 

# Three way sorting - avg distance and flexibility distance
#==========================================================
circuit_metrics.sim.sorted <- sortByTwoIndices.avgDist_and_flex(circuit_metrics.sim) 
write.csv(circuit_metrics.sim.sorted, 
          file = './results/summary.circuits.sim.sortedByAvgDist_flex.csv',  
          row.names = T, quote = F)

