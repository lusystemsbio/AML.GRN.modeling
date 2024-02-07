#!/usr/bin/env Rscript

remove(list = ls()) 

library(sRACIPE)
library(NetAct)
# source('../heatmapSimilarity.updated.R')
# source('../functions.R')
library(gplots)

# Constants 
#----------
circuit_idx <- '0.09-32-0.85'
NO_MODELS <- 10000

# input directory for hS objects obtained from cal.metrics.sim_circuits.R
hS.dir <- '../../networks/circuits.hS/' 

outdir <- './data/'
dir.create(outdir)

# Explore one selected circuit
#-----------------------------
racipe <- readRDS(file = paste('../../networks/circuits.sim/circuit_simulated_', 
                               circuit_idx, '.rds', sep = '') )
names(racipe)
length(names(racipe))


# Plot circuit
#-------------
#sRACIPE::sracipePlotCircuit(racipe, plotToFile = F) 
#sRACIPE:: sracipePlotCircuit(racipe, plotToFile = T) 

class(racipe)
names(racipe)
assays(racipe)

circuit_tpo <- sracipeCircuit(racipe)
length(unique(sort(c(circuit_tpo$Source, circuit_tpo$Target)) ))

# Save the circuit
#----------------- 
outdir <- './data/'
fname.out <- paste(outdir, 'circuit-', circuit_idx, '.tpo' ,sep = '')
write.table(circuit_tpo, file = fname.out, sep = '\t', quote = F, row.names = F)

# extract TFs in the network 
#---------------------------
TFs_in_circuit <- union(circuit_tpo$Source, circuit_tpo$Target)
length(TFs_in_circuit)

fname.out <- paste(outdir, 'TFs-', circuit_idx, '.txt' ,sep = '')
write.table(TFs_in_circuit, file = fname.out, sep = '\t', quote = F, row.names = F)


# Load heatmap similarity object
#-------------------------------
fname.hS <- paste(hS.dir, 'hS_', circuit_idx, '.rds', sep = '')
fname.hS
hS <- readRDS(file = fname.hS)
hS$AvgDist
hS$simulated.cluster.freq
hS$simulated.cluster.freq[2] + hS$simulated.cluster.freq[3]
hS$KL

hS$simulated.cluster.freq
hS$simulated.cluster.freq * NO_MODELS
sum(hS$simulated.cluster.freq * NO_MODELS)
# Save simulation data with cluster information
dim(hS$dataSimulation)
data.sim <- t(hS$dataSimulation)
dim(data.sim) 
cluster.labels <- c(rep('1', hS$simulated.cluster.freq[2]*NO_MODELS),  
                    rep('2', hS$simulated.cluster.freq[3]*NO_MODELS), 
                    rep('3', (as.numeric(hS$simulated.cluster.freq[1]*NO_MODELS))+1)) 

length(cluster.labels) 
sum(cluster.labels=='3') # 642 - Hybrid
sum(cluster.labels=='1') # 2443 - AML
sum(cluster.labels=='2') # 6915 - Untreated

data.sim.labeled <- cbind(cluster.labels, data.sim) 
colnames(data.sim.labeled) <- c('CLUSTER_NO', colnames(data.sim))
fname.out <- paste(outdir, 'racipe.models.wt.labeled.csv', sep = '')
write.csv(data.sim.labeled, fname.out, row.names = FALSE, quote = F)



