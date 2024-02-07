#!/usr/bin/env Rscript
remove(list = ls())

circuit_metrics.sim <- read.csv(file = './results/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)

outdir <- './results/'

# ckt 1
circuit_metrics.sim$Accuracy
prACC = ((sum(circuit_metrics.sim$Accuracy<=circuit_metrics.sim$Accuracy[1]))/nrow(circuit_metrics.sim)*100)
100-prACC # 2.44 ==> 3
prFlex = ((sum(circuit_metrics.sim$flexibility <=circuit_metrics.sim$flexibility[1]))/nrow(circuit_metrics.sim))*100
100-prFlex # 7.52 ==> 8
prACC+prFlex # 10.34 vs 9.3 

prAvg <- (prACC+prFlex)/2
prAvg # 5.17 vs 9.3

trunc(rank(circuit_metrics.sim$Accuracy))/length(circuit_metrics.sim$Accuracy)

head(rank(circuit_metrics.sim$Accuracy))
head(trunc(rank(circuit_metrics.sim$Accuracy)))

tail(rank(circuit_metrics.sim$Accuracy))
tail(trunc(rank(circuit_metrics.sim$Accuracy)))

perc.rank <- function(x) trunc(rank(x))/length(x)

my.df <- as.data.frame(circuit_metrics.sim$Accuracy)
colnames(my.df) <- c('acc') 
my.df <- within(my.df, xr <- perc.rank(acc))


perc.rank.2 <- function(x, xo)  length(x[x <= xo])/length(x)*100
perc.rank.2(circuit_metrics.sim$Accuracy, circuit_metrics.sim$Accuracy) 




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



#======================================
circuit_metrics.sim.list <- readRDS(file = paste0('../samples/comb.metrics.SSIZE.10K/results/circuit_metrics.sim.list.rds'))
names(circuit_metrics.sim.list)


perc.rank.3 <- function(x, xo)  (100-length(x[x <= xo])/length(x)*100)


circuit_m.sim <- circuit_metrics.sim.list$S1

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

circuit_metrics.sim.sorted <- list()



for(sample_name in names(circuit_metrics.sim.list)){ 
  print(sample_name)
  circuit_m.sim <- circuit_metrics.sim.list[[sample_name]]
  prank.acc <- sapply(circuit_m.sim$Accuracy, function (x) perc.rank.3(circuit_metrics.sim$Accuracy, x))
  prank.flex <- sapply(circuit_m.sim$flexibility, function (x) perc.rank.3(circuit_metrics.sim$flexibility, x))

  prank.samples.acc[sample_name] <- prank.acc
  prank.samples.flex[sample_name] <- prank.flex
  
  #break()
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







