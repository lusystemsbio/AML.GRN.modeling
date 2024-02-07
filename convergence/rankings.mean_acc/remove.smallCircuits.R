#!/usr/bin/env Rscript
remove(list = ls()) 

MINIMUM.CIRCUIT.SIZE <- 15

outdir <- './results/'
dir.create(outdir)

circuit_metrics.sim <- read.csv(file = '../../networks/results/summary.circuits.sim.sigVal.csv', row.names = 1)
dim(circuit_metrics.sim) # 555
circuit_metrics.sim <- circuit_metrics.sim[circuit_metrics.sim$Nodes>=MINIMUM.CIRCUIT.SIZE, ]
dim(circuit_metrics.sim) # 532

write.csv(circuit_metrics.sim, file = paste0(outdir, "./summary.circuits.sim.sigVal.csv"), row.names = T, quote = F)
