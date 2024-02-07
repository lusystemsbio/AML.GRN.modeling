#!/usr/bin/env Rscript

remove(list = ls())  

source('../../../networks/functions.R')

datadir <- './results/'

circuit_metrics.sim  <- readRDS(file = paste0(datadir,  './circuit_metrics.sim.sorted.rds'))
names(circuit_metrics.sim)

NO_SAMPLES <- length(names(circuit_metrics.sim))

circuit_metrics <- circuit_metrics.sim$S1

metric.comb <- as.data.frame(matrix(nrow = nrow(circuit_metrics), 
                                    ncol = (length(names(circuit_metrics.sim)))))
#colnames(metric.comb) <- c('SimIdx', names(circuit_metrics.sim))
#colnames(metric.comb) <- c(names(circuit_metrics.sim), 'mean', 'sd')
colnames(metric.comb) <- c(names(circuit_metrics.sim))

circuit_metrics.ordered <- circuit_metrics[order(circuit_metrics$SimIdx), ] 
#metric.comb$SimIdx <- circuit_metrics.ordered$SimIdx 
rownames(metric.comb) <- rownames(circuit_metrics.ordered) 
metric.comb$S1 <- circuit_metrics.ordered$idxBoth

for(sample_name in names(circuit_metrics.sim)[2:NO_SAMPLES]){  
   print(sample_name)
   circuit_metrics <- circuit_metrics.sim[[sample_name]] 
   circuit_metrics.ordered <- circuit_metrics[order(circuit_metrics$SimIdx), ] 
   
   metric.comb[,sample_name] <- circuit_metrics.ordered$idxBoth
}

class(metric.comb$S1)

avg.samples <- apply(metric.comb, 1, mean)
sd.samples <- apply(metric.comb, 1, sd)
length(avg.samples)
length(sd.samples)

metric.comb$mean <- avg.samples
metric.comb$sd <- sd.samples

metric.comb.sorted <- metric.comb[order(metric.comb$mean),]

cat(head(rownames(metric.comb.sorted)))

write.csv(metric.comb.sorted, file = paste0(datadir, 'metric.comb.sorted.csv'), 
          row.names = T)
