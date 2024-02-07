#!/usr/bin/env Rscript

remove(list = ls()) 

# Load circuit summary 
#-------------------- 
#datadir <- '../rankingsBypercentile/results/'
#datadir <- './results/'
#circuit_metrics.sim <- read.csv(file = paste0(datadir , "./summary.circuits.sim.sortedByAcc_flex.csv"), row.names = 1)
circuit_metrics.sim <- read.csv(file = "./results/summary.circuits.sim.sigVal.csv", row.names = 1)

# Load mean accuracy
#--------------------  
mean_sd.df <- read.csv(file = '../accuracy.fluctuation/para.all/results/accuracy.mean_sd.csv', row.names = 1)

# bring the circutis in the same order 
sum(rownames(circuit_metrics.sim)==rownames(mean_sd.df))
circuit_metrics.sim <- circuit_metrics.sim[rownames(mean_sd.df),]
sum(rownames(circuit_metrics.sim)==rownames(mean_sd.df))

# plot(circuit_metrics.sim$Accuracy)
# points(mean_sd.df$mean, col='red')
# 
# plot(circuit_metrics.sim$Accuracy[1:50], ylim = c(0, 1))
# points(mean_sd.df$mean[1:50], col='red')
# 
# plot(circuit_metrics.sim$Accuracy[1:10], ylim = c(0, 1))
# points(mean_sd.df$mean[1:10], col='red')

# replace accuracy with mean accuracy
circuit_metrics.sim$Accuracy <- mean_sd.df$mean

# add sd of accuracies
library("tibble")
circuit_metrics.sim  <- add_column(circuit_metrics.sim, Accuracy_sd=mean_sd.df$sd, .after = 'Accuracy')

write.csv(format(circuit_metrics.sim, digits=3), file = paste0("./results/summary.circuits.sim.sigVal.csv"), row.names = T, quote = F) 
