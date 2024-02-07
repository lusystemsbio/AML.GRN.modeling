#!/usr/bin/env Rscript
remove(list = ls())

outdir <- './results/'
circuit_metrics.sim <- read.csv(file = './results/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)

circuit_metrics.sim.list <- readRDS(file = paste0('../samples/comb.metrics.SSIZE.10K/results/circuit_metrics.sim.list.rds'))
names(circuit_metrics.sim.list)

perc.rank.3 <- function(x, xo)  (100-length(x[x <= xo])/length(x)*100)
#perc.rank.3 <- function(x, xo)  length(x[x >= xo])/length(x)*100

# create data structures
circuit_m.sim <- circuit_metrics.sim.list$S1

# check whetehr they are the same circuits
topCkts <- rownames(circuit_metrics.sim)[1:10]
sum(rownames(circuit_m.sim) == topCkts)

prank.samples.acc <- as.data.frame(matrix(nrow = nrow(circuit_m.sim), 
                                          ncol = (length(names(circuit_metrics.sim.list))+1)))
rownames(prank.samples.acc) <- row.names(circuit_metrics.sim)[1:10]
colnames(prank.samples.acc) <- c('nodes', names(circuit_metrics.sim.list))
prank.samples.acc$nodes <- circuit_m.sim$Nodes

prank.samples.flex <- as.data.frame(matrix(nrow = nrow(circuit_m.sim), 
                                          ncol = (length(names(circuit_metrics.sim.list))+1)))
rownames(prank.samples.flex) <- row.names(circuit_metrics.sim)[1:10]
colnames(prank.samples.flex) <- c('nodes', names(circuit_metrics.sim.list))
prank.samples.flex$nodes <- circuit_m.sim$Nodes

# perform ranking for each sample
for(sample_name in names(circuit_metrics.sim.list)){ 
  print(sample_name)
  circuit_m.sim <- circuit_metrics.sim.list[[sample_name]]
  # prank.acc <- sapply(circuit_m.sim$Accuracy, function (x) perc.rank.3(circuit_metrics.sim$Accuracy, x))
  # prank.flex <- sapply(circuit_m.sim$flexibility, function (x) perc.rank.3(circuit_metrics.sim$flexibility, x))

  prank.samples.acc[sample_name] <- sapply(circuit_m.sim$Accuracy, function (x) perc.rank.3(circuit_metrics.sim$Accuracy, x))
  prank.samples.flex[sample_name] <- sapply(circuit_m.sim$flexibility, function (x) perc.rank.3(circuit_metrics.sim$flexibility, x))
  #print(sum(rownames(circuit_m.sim) == topCkts))
}

prank.samples.acc$mean <- apply(prank.samples.acc[2:ncol(prank.samples.acc)], 1, mean)
prank.samples.acc$sd <- apply(prank.samples.acc[2:ncol(prank.samples.acc)], 1, sd)

prank.samples.flex$mean <- apply(prank.samples.flex[2:ncol(prank.samples.flex)], 1, mean)
prank.samples.flex$sd <- apply(prank.samples.flex[2:ncol(prank.samples.flex)], 1, sd)


prank.acc.bg <- sapply(circuit_metrics.sim$Accuracy[1:10], function (x) perc.rank.3(circuit_metrics.sim$Accuracy, x))
prank.flex.bg <- sapply(circuit_metrics.sim$flexibility[1:10], function (x) perc.rank.3(circuit_metrics.sim$flexibility, x))


prank.samples.acc$bg <- prank.acc.bg
prank.samples.flex$bg <- prank.flex.bg


prank.samples.avg <- cbind(prank.samples.flex$nodes, 
                           (prank.samples.acc$mean + prank.samples.flex$mean)/2, 
                           (prank.samples.acc$bg+prank.samples.flex$bg)/2) 

rownames(prank.samples.avg) <- rownames(prank.samples.flex)
colnames(prank.samples.avg) <- c('nodes', 'avg.rank.samples', 'avg.rank.bg')

write.csv(format(prank.samples.acc, digits = 3), file = './results/prank.samples.acc.csv', row.names = T, quote = F) 
write.csv(format(prank.samples.flex, digits = 3), file = './results/prank.samples.flex.csv', row.names = T, quote = F) 
write.csv(format(prank.samples.avg, digits = 3), file = './results/prank.samples.avg.csv', row.names = T, quote = F) 

