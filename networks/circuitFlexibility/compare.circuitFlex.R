#!/usr/bin/env Rscript

remove(list = ls()) 
circuitFlex.euc <- read.csv(file = './data/circuit.flexibility.csv')
circuitFlex.bc <- read.csv(file = './data/flexibility_BC.csv')

circuitFlex.euc <- circuitFlex.euc[order(circuitFlex.euc$circuit.idx, decreasing = TRUE), ]
circuitFlex.bc <- circuitFlex.bc[order(circuitFlex.bc$circuit.idx, decreasing = TRUE), ]
rownames(circuitFlex.bc) <- circuitFlex.bc$circuit.idx

cor(circuitFlex.euc$flexibility, circuitFlex.bc$avg.dist.kd)
cor(circuitFlex.euc$flexibility, circuitFlex.bc$avg.dist.oe)
cor(circuitFlex.bc$avg.dist.kd, circuitFlex.bc$avg.dist.oe) 
cor(circuitFlex.euc$flexibility, circuitFlex.bc$flexibility)

circuit_metrics.sim <- read.csv(file = '../results/flex.euc/summary.circuits.sim.sorted.csv', row.names = 1)

circuit_metrics.new <- circuit_metrics.sim[, c(1:25)]
circuitFlex.bc <- circuitFlex.bc[rownames(circuit_metrics.sim),] 

circuit_metrics.new$flexibility.bc <- circuitFlex.bc$avg.dist.kd

circuit_metrics.new <- circuit_metrics.new[order(circuit_metrics.new$flexibility.bc, decreasing = TRUE), ]
circuit_metrics.new$idxFlexbc <- c(1:dim(circuit_metrics.new)[1])
circuit_metrics.new$idxTrio <- circuit_metrics.new$idxAccuracy+circuit_metrics.new$idxAvgDist + circuit_metrics.new$idxFlexbc

circuit_metrics.new <- circuit_metrics.new[order(circuit_metrics.new$idxTrio, decreasing = FALSE), ]

circuit_metrics.new[1:10, 1:6] 
circuit_metrics.new[219:228, 1:6] 
sort(circuit_metrics.new$Nodes[1:10])
unique(sort(circuit_metrics.new$Nodes[219:228]))

# how many top circuits are found to be common
#----------------------------------------------
topCircuits.euc <- rownames(circuit_metrics.sim[1:10, 1:6])
topCircuits.bc <- rownames(circuit_metrics.new[1:10, 1:6])

topCircuits.common <- intersect(topCircuits.euc, topCircuits.bc)
length(topCircuits.common)

circuit_metrics.common <- circuit_metrics.new[topCircuits.common, ]
circuit_metrics.common <- circuit_metrics.sim[topCircuits.common, ]
circuit_metrics.common[1:6, 1:6] 
