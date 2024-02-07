#!/usr/bin/env Rscript

remove(list = ls())

# Parameters
#-----------
outdir <- './degs/'
dir.create(outdir)

library(NetAct)

# Load DE results
#----------------
fname.de.results <- '../../data.tfs/de.results.rda'
load(fname.de.results) # loads object de.results
names(de.results)

# GSEA using only the network genes
#-----------------------------------
#circuit_metrics.sim <- read.csv(file = '../networks/results/summary.circuits.sim.sorted.csv', row.names = 1)
#circuit_idx <- rownames(circuit_metrics.sim)[1]

circuit_idx <- "0.09-32-0.85"
circuit_idx


fname.hS <- paste('../../networks/circuits.hS/hS_', circuit_idx, '.rds', sep = '')
hS <- readRDS(file = fname.hS)
class(hS)
names(hS)
data.sim <- hS$dataSimulation
class(data.sim)
data.ref <- hS$dataReference
dim(data.ref)

TFs_in_circuit <- rownames(data.ref)
length(TFs_in_circuit) 

# Obtain target DB relevant to this circuit
#------------------------------------------
regDB.list <- readRDS('../../databases/targetDB.list.rds') 
regDB.feature.ratio <- feature_ratio_cutoff <- strsplit(circuit_idx, '-', 2)[[1]][1]  
regDB.comb <- regDB.list[[regDB.feature.ratio]] 
length(names(regDB.comb)) # 455 

regDB.circuit <- lapply(TFs_in_circuit, function(tf) regDB.comb[[tf]]) 
length(regDB.circuit) 
names(regDB.circuit) <- TFs_in_circuit
head(regDB.circuit)

# combine the TFs and their targets: 
TFsPlusTargets_in_cicuit <- unique(c(names(regDB.circuit), unlist(regDB.circuit))) 
length(TFsPlusTargets_in_cicuit) # 3319

# Common genes between TFsPlusTargets_in_cicuit and de.results 
#-------------------------------------------------------------
genes.common <- intersect(TFsPlusTargets_in_cicuit, de.results$degs)
length(genes.common) # 1468 

fname.out <- paste(outdir, 'DEGs_in_circuit.txt', sep = '')
write.table(genes.common, file = fname.out, 
            col.names = F, row.names = F, quote = F)

