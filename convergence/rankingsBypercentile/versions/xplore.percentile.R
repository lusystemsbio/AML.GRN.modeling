#!/usr/bin/env Rscript
remove(list = ls())

circuit_metrics.sim <- read.csv(file = './results/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)


# ckt 1
circuit_metrics.sim$idxAccuracy[1]
prACC = ((sum(circuit_metrics.sim$idxAccuracy<=circuit_metrics.sim$idxAccuracy[1]))/nrow(circuit_metrics.sim)*100)
prACC # 2.63
prFlex = ((sum(circuit_metrics.sim$idxFlexibility<=circuit_metrics.sim$idxFlexibility[1]))/nrow(circuit_metrics.sim))*100
prFlex # 7.71
prACC+prFlex # 10.34 vs 9.3 

prAvg <- (prACC+prFlex)/2
prAvg # 5.17 vs 9.3




# ckt 2
circuit_metrics.sim$idxAccuracy[2]
prACC = ((sum(circuit_metrics.sim$idxAccuracy<=circuit_metrics.sim$idxAccuracy[2]))/nrow(circuit_metrics.sim)*100)
prACC
prFlex = ((sum(circuit_metrics.sim$idxFlexibility<=circuit_metrics.sim$idxFlexibility[2]))/nrow(circuit_metrics.sim))*100
prFlex
prACC+prFlex  # 17.29323 vs 11.7 

prAvg <- (prACC+prFlex)/2
prAvg # 8.65 vs 11.7
 

# ckt 3
circuit_metrics.sim$idxAccuracy[3]
prACC = ((sum(circuit_metrics.sim$idxAccuracy<=circuit_metrics.sim$idxAccuracy[3]))/nrow(circuit_metrics.sim)*100)
prACC
prFlex = ((sum(circuit_metrics.sim$idxFlexibility<=circuit_metrics.sim$idxFlexibility[3]))/nrow(circuit_metrics.sim))*100
prFlex
prACC+prFlex # 19.55 vs 8.1
prAvg <- (prACC+prFlex)/2
prAvg # 9.77  vs 8.1




circuits.small <- circuit_metrics.sim[circuit_metrics.sim$Nodes<=20, ]



# samples 10K
circuit_metrics.sim.samples <- readRDS(file = paste0('../samples/comb.metrics.SSIZE.10K/rankedCircuits.acc.flex/circuit_metrics.sim.sorted.rds'))
s1 <- circuit_metrics.sim.samples$S1
s2 <- circuit_metrics.sim.samples$S2
s3 <- circuit_metrics.sim.samples$S3
s4 <- circuit_metrics.sim.samples$S4
s5 <- circuit_metrics.sim.samples$S5
s6 <- circuit_metrics.sim.samples$S6
s7 <- circuit_metrics.sim.samples$S7
s8 <- circuit_metrics.sim.samples$S8
s9 <- circuit_metrics.sim.samples$S9
s10 <- circuit_metrics.sim.samples$S10





# mean and sds 
metric.comb.sorted <- read.csv(file = '../samples/comb.metrics.SSIZE.10K/rankedCircuits.acc.flex/metric.comb.sorted.csv', row.names = 1)
