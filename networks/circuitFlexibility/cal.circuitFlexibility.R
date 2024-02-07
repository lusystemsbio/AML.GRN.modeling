#!/usr/bin/env Rscript

remove(list = ls())

PERCENT_REDUCTION <- 90

library(sRACIPE)
source('./functions.flexibility.R')

# Load circuits
circuits <- readRDS('../circuits/circuits.uniq.rds')

names(circuits)

fname.list <- list() 
for(idx in names(circuits)){
   #print(idx)
   fname.list[[idx]] <- paste('./circuits.sim/circuit_simulated_', idx, '.rds', sep = '')
}
names(fname.list) <- names(circuits)
 
fname.list[1]

flexibility.df <- as.data.frame(matrix(nrow = length(circuits), ncol = 2)) 
colnames(flexibility.df) <- c('circuit.idx', 'flexibility')
rownames(flexibility.df) <- names(circuits)

for(circuit_idx in names(circuits)){
   #circuit_idx <- names(circuits)[1]
   #circuit_idx <- '0.16-16-0.8' #'0.07-20-0.8' #'0.07-12-0.8' #'0.07-12-0.85' #'0.06-20-0.8' #'0.06-16-0.8' #'0.06-12-0.8'
   print(circuit_idx)
   
   # Obtain simulated data:
   racipe <- readRDS(file = paste('../circuits.sim/circuit_simulated_', 
                                  circuit_idx, '.rds', sep = '')) 
   racipe <- sRACIPE::sracipeNormalize(racipe) 
   
   # KD by subsetting:
   racipe.kd <- sRACIPE::sracipeKnockDown(racipe, plotToFile = FALSE,
                                          plotBarPlot = FALSE, #TRUE, 
                                          plotHeatmap = FALSE, 
                                          reduceProduction = 10 #3 #2 #(100-PERCENT_REDUCTION)
                                          )
   avg.dist <- calDistance(racipe.kd)
   flexibility.df[circuit_idx, ] <- c(circuit_idx, avg.dist)
   #break()
}

#sRACIPE::sracipePlotCircuit(racipe, plotToFile = FALSE)

fname <- './circuit.flexibility.csv'
write.csv(flexibility.df, file = fname, quote = FALSE, row.names = FALSE)

