#!/usr/bin/env Rscript

remove(list = ls()) 

#source('../../../networks/functions.R')
source('./functions.R')

datadir <- './results/'
outdir <- './rankedCircuits.acc.flex/'
dir.create(outdir)

circuit_metrics.sim.list <- readRDS(file = paste0(datadir,  './circuit_metrics.sim.list.rds'))
names(circuit_metrics.sim.list)

circuit_metrics.sim.sorted <- list()
# sample_name <- names(circuit_metrics.sim.list)[1]
# sample_name <- names(circuit_metrics.sim.list)[2]
# sample_name <- names(circuit_metrics.sim.list)[3] 
# sample_name <- names(circuit_metrics.sim.list)[4] 
# sample_name <- names(circuit_metrics.sim.list)[5] 
# sample_name <- names(circuit_metrics.sim.list)[6] 
# sample_name <- names(circuit_metrics.sim.list)[7] 
# sample_name <- names(circuit_metrics.sim.list)[8] 
# sample_name <- names(circuit_metrics.sim.list)[9] 
# sample_name <- names(circuit_metrics.sim.list)[10] 

for(sample_name in names(circuit_metrics.sim.list)){
   print(sample_name)
   circuit_metrics.sim <- circuit_metrics.sim.list[[sample_name]]
   #circuit_metrics.sim.s <- sortByTwoIndices(circuit_metrics.sim)
   #circuit_metrics.sim.sorted[[sample_name]]  <- sortByTwoIndices(circuit_metrics.sim.list[[sample_name]])
   
   circuit_metrics.sim.s <- sortByTwoIndices.acc_and_flex(circuit_metrics.sim)
   circuit_metrics.sim.sorted[[sample_name]]  <- sortByTwoIndices.acc_and_flex(circuit_metrics.sim.list[[sample_name]]) 
   
   # break()
}
names(circuit_metrics.sim.sorted)

saveRDS(circuit_metrics.sim.sorted, file = paste0(outdir,  './circuit_metrics.sim.sorted.rds'))
