#!/usr/bin/env Rscript

remove(list = ls())

PERCENT_REDUCTION <- 90

library(sRACIPE)
source('./sracipeGeneImp.R')

# Load circuits
circuits <- readRDS('../circuits/circuits.uniq.rds')

names(circuits)

flexibility.df <- as.data.frame(matrix(nrow = length(circuits), ncol = 3)) 
colnames(flexibility.df) <- c('circuit.idx', 'avg.dist.kd', 'avg.dist.oe')
rownames(flexibility.df) <- names(circuits)

for(circuit_idx in names(circuits)){
   #circuit_idx <- names(circuits)[1]
   print(circuit_idx)

   # Obtain simulated data:
   racipe <- readRDS(file = paste('../circuits.sim/circuit_simulated_', 
                                  circuit_idx, '.rds', sep = '')) 
   racipe <- sRACIPE::sracipeNormalize(racipe) 
   
   # KD by subsetting:
   distByGene <- sracipeGeneImp(racipe = racipe, kdPercent=10)
   flexibility.df[circuit_idx, ] <- c(circuit_idx, as.numeric(colMeans(distByGene)))
   #break()
}

#sRACIPE::sracipePlotCircuit(racipe, plotToFile = FALSE)

# calculate average of the two distances as the final circuit flexibility:
circuitFlex.bc$flexibility <- (circuitFlex.bc$avg.dist.kd+circuitFlex.bc$avg.dist.oe)/2

fname <- './data/flexibility_BC.csv'
write.csv(flexibility.df, file = fname, quote = FALSE, row.names = FALSE)
