#!/usr/bin/env Rscript

remove(list = ls())  

datadir <- './rankedCircuits.acc.flex/'

circuit_metrics.sim  <- readRDS(file = paste0(datadir,  './circuit_metrics.sim.sorted.rds'))
names(circuit_metrics.sim)

NO_SAMPLES <- length(names(circuit_metrics.sim))

circuit_metrics <- circuit_metrics.sim$S1

metric.comb <- as.data.frame(matrix(nrow = nrow(circuit_metrics), 
                                    ncol = (length(names(circuit_metrics.sim))+1)))
colnames(metric.comb) <- c('Nodes' ,names(circuit_metrics.sim))

circuit_metrics.ordered <- circuit_metrics[order(circuit_metrics$SimIdx), ] 
rownames(metric.comb) <- rownames(circuit_metrics.ordered)
metric.comb$Nodes <- circuit_metrics.ordered$Nodes
metric.comb$S1 <- circuit_metrics.ordered$idxFlexibility

for(sample_name in names(circuit_metrics.sim)[2:NO_SAMPLES]){  
   print(sample_name)
   circuit_metrics <- circuit_metrics.sim[[sample_name]] 
   circuit_metrics.ordered <- circuit_metrics[order(circuit_metrics$SimIdx), ] 
   metric.comb[,sample_name] <- circuit_metrics.ordered$idxFlexibility
}

class(metric.comb$S1)

# avg.samples <- apply(metric.comb, 1, mean)
# sd.samples <- apply(metric.comb, 1, sd)
avg.samples <- apply(metric.comb[2:ncol(metric.comb)], 1, mean)
sd.samples <- apply(metric.comb[2:ncol(metric.comb)], 1, sd)
length(avg.samples)
length(sd.samples)

metric.comb$mean <- avg.samples
metric.comb$sd <- sd.samples

metric.comb.sorted <- metric.comb[order(metric.comb$mean),]

cat(head(rownames(metric.comb.sorted)))

write.csv(format(metric.comb.sorted, digits=3), file = paste0(datadir, 'metric.idxFlex.sorted.csv'), 
          row.names = T, quote = F)
