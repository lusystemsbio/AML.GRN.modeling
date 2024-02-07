#!/usr/bin/env Rscript

remove(list = ls()) 

library(sRACIPE)

# Constants 
#------------
outdir <- './tfs.annotated/'

# Obtain circuit index for the top circuit
#-----------------------------------------
circuit_metrics <- read.csv(file = '../networks/results/summary.circuits.sim.sorted.csv', row.names = 1)
circuit_idx <- rownames(circuit_metrics)[1]
circuit_idx # "0.07-32-0.75" - top circuit in all topologies

# Load racipe object 
#------------------------------
racipe <- readRDS(file = paste('../networks/circuits.sim/circuit_simulated_', 
                               circuit_idx, '.rds', sep = '') )
names(racipe)
length(names(racipe))

# Extract circt
#--------------
circuit_tpo <- sracipeCircuit(racipe)
length(unique(sort(c(circuit_tpo$Source, circuit_tpo$Target)) ))

# Save the circuit
#-----------------
fname.out <- paste(outdir, 'circuit-', circuit_idx, '.tpo' ,sep = '')
write.table(circuit_tpo, file = fname.out, sep = '\t', quote = F, row.names = F)
