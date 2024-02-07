#!/usr/bin/env Rscript

remove(list = ls())


library(sRACIPE)

# top circuits: 1 ~ 10
# bottom circuits: 219 ~ 228

circuit_metrics <- read.csv(file = '../results/summary.circuits.sim.sorted.csv', row.names = 1)
circuit_idx <- rownames(circuit_metrics)[1] 
circuit_idx  

racipe <- readRDS(file = paste('../circuits.sim/circuit_simulated_', 
                               circuit_idx, '.rds', sep = '') )


#sRACIPE::sracipePlotCircuit(racipe, plotToFile = FALSE)

# Find the common nodes between the best and the worst circuit 
#-------------------------------------------------------------
circuit.tpo <- sRACIPE::sracipeCircuit(racipe)
nodes1 <- unique(sort(union(circuit.tpo$Source, circuit.tpo$Target)))
nodes228 <- unique(sort(union(circuit.tpo$Source, circuit.tpo$Target)))

length(nodes1)
length(nodes228)

nodes.common <- intersect(nodes1, nodes228)
length(nodes.common)
cat(nodes.common)
