#!/usr/bin/env Rscript
remove(list = ls()) 

circuits.u <- readRDS('./circuits/circuits.uniq.rds')
fnames <- paste('circuit_simulated_', names(circuits.u), '.rds', sep = '')
head(fnames)
length(fnames)
filepaths <- paste('./circuits.sim.both/', fnames, sep = '')
filepaths
#file.copy(from = filepaths, to='./circuits.sim/')
