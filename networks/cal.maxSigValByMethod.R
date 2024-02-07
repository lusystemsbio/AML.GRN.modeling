#!/usr/bin/env Rscript

remove(list = ls()) 

# PARAMETERS

# NUM_CLUSTERS: user supplied information about 
# how many clusters in the data
NUM_CLUSTERS <- 2

library(sRACIPE)
source('./functions.R')
library(dplyr)

outdir <- './results/'
dir.create(outdir)

# Load circuit summary
#-------------------- 
circuit_metrics.sim <- read.csv(file = './results/summary.circuits.sim.flex.csv', row.names = 1)

# Load circuits
circuits <- readRDS('./circuits/circuits.uniq.rds')
names(circuits)
length(names(circuits))
fname.list <- list() 
for(idx in names(circuits)){
   print(idx)
   fname.list[[idx]] <- paste('./circuits.sim/circuit_simulated_', idx, '.rds', sep = '')
}
names(fname.list) <- names(circuits)


# select TOP n TFs from each of the THREE methods: NetAct, MARINa, and RI
#------------------------------------------------------------------------
coreTFs.list <- readRDS('../tfSets/data/coreTFs.rds')
names(coreTFs.list)

# Construct possible circuit set for the TFs based on each TF-target DB in the TF-target DB set
#---------------------------------------------------------------------------------------------- 
targetDB.list <- readRDS('../databases/targetDB.list.rds')
names(targetDB.list)

# Calculate heatmap similarity between data.REF and data.sim 
#------------------------------------------------------------
for(circuit_idx in names(circuits)){
   #circuit_idx <- names(circuits)[1]
   print(circuit_idx)
   #break()
   # obtain TF activities using NetAct:
   fr <- strsplit(circuit_idx, split = '-')[[1]][1]
   top.TFs.count <- strsplit(circuit_idx, split = '-')[[1]][2]  
   
   #coreTFs <- coreTFs.list[[top.TFs.count]] 
   coreTFs <- coreTFs.list[[fr]][[top.TFs.count]][['COMB']]
   targetDB = targetDB.list[[fr]] 
   length(coreTFs)
   coreTFs <- intersect(coreTFs, names(targetDB))
   length(coreTFs) 

   # obtain simulated data:
   racipe <- readRDS(file = paste('./circuits.sim/circuit_simulated_', 
                                        circuit_idx, '.rds', sep = '')) 
   racipe <- sracipeNormalize(racipe)  

   # obtain the significance of TFs in the circuit: 
   circuit.df <- sRACIPE::sracipeCircuit(racipe)
   circuit.nodes <- unique(sort(union(circuit.df[,1], circuit.df[,2]))) 
   coreTFs.Netact <- coreTFs.list[[fr]][[top.TFs.count]][['Netact']]
   coreTFs.MARINa <- coreTFs.list[[fr]][[top.TFs.count]][['MARINa']]
   coreTFs.RI <- coreTFs.list[[fr]][[top.TFs.count]][['RI']]  

   maxSigval.Netact <- obtain.sigval.Netact_or_MARINa(coreTFs.Method=coreTFs.Netact, 
                                                      circuit.nodes=circuit.nodes) 
   maxSigval.MARINa <- obtain.sigval.Netact_or_MARINa(coreTFs.Method=coreTFs.MARINa, 
                                                      circuit.nodes=circuit.nodes) 
   minSigval.RI <- obtain.sigval.RI(coreTFs.Method=coreTFs.RI, circuit.nodes=circuit.nodes) 
   circuit_metrics.sim[circuit_idx, "maxSigval.Netact"] <- maxSigval.Netact
   circuit_metrics.sim[circuit_idx, "maxSigval.MARINa"] <- maxSigval.MARINa
   circuit_metrics.sim[circuit_idx, "minSigval.RI"] <- minSigval.RI
   #break()
}

write.csv(circuit_metrics.sim, file = paste0(outdir, "./summary.circuits.sim.sigVal.csv"),
          row.names = T, quote = F) 
