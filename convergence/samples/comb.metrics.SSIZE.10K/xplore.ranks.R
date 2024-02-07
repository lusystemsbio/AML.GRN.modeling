#!/usr/bin/env Rscript
remove(list = ls()) 

circuit_metrics.sim.sorted <- readRDS(file = paste0('./rankedCircuits.acc.flex/circuit_metrics.sim.sorted.rds'))

names(circuit_metrics.sim.sorted)

s1 <- circuit_metrics.sim.sorted$S1
s2 <- circuit_metrics.sim.sorted$S2
s3 <- circuit_metrics.sim.sorted$S3
s4 <- circuit_metrics.sim.sorted$S4
s5 <- circuit_metrics.sim.sorted$S5
s6 <- circuit_metrics.sim.sorted$S6
s7 <- circuit_metrics.sim.sorted$S7
s8 <- circuit_metrics.sim.sorted$S8
s9 <- circuit_metrics.sim.sorted$S9
s10 <- circuit_metrics.sim.sorted$S10

circuit_metrics.sim.list <- readRDS(file = paste0('./results/circuit_metrics.sim.list.rds'))
names(circuit_metrics.sim.list)

circuit_metrics.sim.list$S1

s1 <- circuit_metrics.sim.list$S1
