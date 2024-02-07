#!/usr/bin/env Rscript

remove(list = ls())
NO_TOP_CIRCUITS <- 10

library(sRACIPE)

# Load top circuit ids 
#---------------------
circuit_summary <- read.csv(file = './results/summary.circuits.sim.sorted.csv', row.names = 1)
circuit_idx <- rownames(circuit_summary) #[1]
topCircuits <- circuit_idx[1:NO_TOP_CIRCUITS]
#rm(circuit_summary)
#circuit_idx


# Load geneSets
#--------------
coreTFs.list <- readRDS('../tfSets/data/coreTFs.rds')
names(coreTFs.list)

topCircuits

for(fr in names(topCircuits)){
   coreTFs.byTOPtfCount <- coreTFs.list[[fr]] 
   for(top.TFs.count in names(coreTFs.byTOPtfCount)){
      print(top.TFs.count) 
      coreTFs.byMethod <- coreTFs.byTOPtfCount[[top.TFs.count]] 
      
      break()
   }
   break() 
} 


tfCount.byMetod <- as.data.frame(matrix(nrow = dim(circuit_summary)[1], ncol = 3))
rownames(tfCount.byMetod) <- rownames(circuit_summary) 
colnames(tfCount.byMetod) <- c('tfCount.Netact', 'tfCount.MARINa', 'tfCount.RI')
for(circuit_idx in rownames(circuit_summary)){
   print(circuit_idx) 
   # Obtain simulated data:
   racipe <- readRDS(file = paste('./circuits.sim/circuit_simulated_', 
                                  circuit_idx, '.rds', sep = ''))  
   
   circuit.tpo <- sRACIPE::sracipeCircuit(racipe) 
   tfs_on_circuit <- unique(sort(union(circuit.tpo$Source, circuit.tpo$Target))) 
   fr <- strsplit(circuit_idx, split = '-', fixed = TRUE)[[1]][1] 
   top.TFs.count <- strsplit(circuit_idx, split = '-', fixed = TRUE)[[1]][2] 
   coreTFs.byMethod <- coreTFs.list[[fr]][[top.TFs.count]] 
   tfCount.byMetod[circuit_idx, ] <- c(length(intersect(tfs_on_circuit, coreTFs.byMethod[["Netact"]])), 
                                       length(intersect(tfs_on_circuit, coreTFs.byMethod[["MARINa"]])),
                                       length(intersect(tfs_on_circuit, coreTFs.byMethod[["RI"]])))
   
   #break()
}

tfCount.byMetod

tmp.df <- tfCount.byMetod 
tfCount.byMetod <- cbind(circuit_summary[1:6], tmp.df)
colnames(tfCount.byMetod) <- c(colnames(circuit_summary[1:6]), colnames(tmp.df))

outdir <- './results/'
fname <- paste(outdir, 'tfCount.byMetod.csv', sep = '')  
write.csv(tfCount.byMetod, file = fname, quote = FALSE, row.names = FALSE)

