#!/usr/bin/env Rscript

remove(list = ls())

PERCENT_REDUCTION <- 90

library(sRACIPE)
source('./sracipeGeneImp.R')


# Load circuits
circuits <- readRDS('../circuits/circuits.uniq.rds')
circuit_idx <- '0.16-28-0.85'
sum(names(circuits)==circuit_idx)

# Obtain simulated data:
racipe <- readRDS(file = paste('../circuits.sim/circuit_simulated_', 
                               circuit_idx, '.rds', sep = '')) 
racipe <- sRACIPE::sracipeNormalize(racipe) 

# KD by subsetting:
distByGene <- sracipeGeneImp(racipe = racipe, kdPercent=10)



