#!/usr/bin/env Rscript

remove(list = ls()) 


circuit_metrics.sim.new <- read.csv(file = '../../networks/results/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)
circuit_metrics.sim.old <- read.csv(file = '../../networks/results.old/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)
circuit_metrics.sim.old <- circuit_metrics.sim.old[circuit_metrics.sim.old$Nodes>=15, ]

circuit_metrics.sim.old.2 <- circuit_metrics.sim.old[rownames(circuit_metrics.sim.new), ]

sum(rownames(circuit_metrics.sim.old.2) == rownames(circuit_metrics.sim.new))

sum(circuit_metrics.sim.new$flexibility==circuit_metrics.sim.old.2$flexibility)


for(idx in rownames(circuit_metrics.sim.new)[1:10]){
  print(idx) 
  racipe <- readRDS(file = paste('../networks/circuits.sim/circuit_simulated_', idx, '.rds', sep = '') ) 
  
}



