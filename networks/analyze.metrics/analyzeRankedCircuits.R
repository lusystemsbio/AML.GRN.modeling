#!/usr/bin/env Rscript

remove(list = ls())

# Resources: 
# 1. http://r-statistics.co/Linear-Regression.html 

circuit_metrics.sim <- read.csv(file = '../results/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)

max(circuit_metrics.sim$Nodes)
library(sRACIPE)
library(NetAct)
source('../heatmapSimilarity.updated.R')
source('../functions.R')
library(gplots)


# Constants 
#------------
CEX.MAIN <- 1.5
CEX.LAB <- 1.5
CEX.AXIS <- 1.5


# input directory for hS objects obtained from cal.metrics.sim_circuits.R
hS.dir <- '../circuits.hS/' 

outdir <- './figs/'
dir.create(outdir)


# Correlation values between Accuracy, AvgDist and KL distance in the sorted data
#--------------------------------------------------------------------------------
cor(circuit_metrics.sim$AvgDist, circuit_metrics.sim$Accuracy) 
cor(circuit_metrics.sim$AvgDist, circuit_metrics.sim$KLdist)
cor(circuit_metrics.sim$Accuracy, circuit_metrics.sim$KLdist)


cor(circuit_metrics.sim$Accuracy, circuit_metrics.sim$flexibility)
cor(circuit_metrics.sim$AvgDist, circuit_metrics.sim$flexibility)
cor(circuit_metrics.sim$KLdist, circuit_metrics.sim$flexibility)

cor(circuit_metrics.sim$TopTFs, circuit_metrics.sim$Nodes)

# Comparing Accuracy, AvgDist and KL distance in the sorted data
#---------------------------------------------------------------
fname_fig <- paste0(outdir, 'average_distance.after_sorting.pdf')
pdf(file=fname_fig, width=6, height=8, paper = "special", onefile = TRUE)
par(mfrow=c(3,1))
par(mar=c(2.0, 5.5, 1.5, 2)) # bottom, left, top, right 
plot(circuit_metrics.sim$idxBoth, circuit_metrics.sim$Accuracy, xaxt="n", ylab='', cex.axis=CEX.AXIS)
title(ylab='Accuracy', cex.lab=CEX.LAB)
plot(circuit_metrics.sim$idxBoth, circuit_metrics.sim$AvgDist, xaxt="n", ylab='', cex.axis=CEX.AXIS)
title(ylab='Average distance', cex.lab=CEX.LAB)
par(mar=c(5.0, 5.5, 1.5, 2))
plot(circuit_metrics.sim$idxBoth, circuit_metrics.sim$KLdist, xlab='', ylab='', cex.axis=CEX.AXIS)
title(xlab='Sort index (combined)', ylab='KL distance', cex.lab=CEX.LAB)
dev.off()


# Histograms for the normalized metrics 
#--------------------------------------
# normalize THREE metricx
#-------------------------
avg.distance.normed <- (circuit_metrics.sim$AvgDist - min(circuit_metrics.sim$AvgDist))/(max(circuit_metrics.sim$AvgDist)-min(circuit_metrics.sim$AvgDist))
accuracy.normed <- (circuit_metrics.sim$Accuracy - min(circuit_metrics.sim$Accuracy))/(max(circuit_metrics.sim$Accuracy)-min(circuit_metrics.sim$Accuracy))
KLdistance.normed <- (circuit_metrics.sim$KLdist-min(circuit_metrics.sim$KLdist))/(max(circuit_metrics.sim$KLdist)-min(circuit_metrics.sim$KLdist))


BREAKS <- 60
par(mfrow=c(3,1))
par(mar=c(2.0, 5.5, 1.5, 2)) # bottom, left, top, right 
h1 <- hist(accuracy.normed, breaks = BREAKS, main = 'Accuracy') 
h2 <- hist(avg.distance.normed, breaks = BREAKS, main = 'Average distance') 
h3 <- hist(KLdistance.normed, breaks = BREAKS, main = 'KL distance') 



XLIMIT <- c(0, 1)
fname_fig <- paste0(outdir, 'histogram-metrics.pdf')
pdf(file=fname_fig, width=6, height=8, paper = "special", onefile = TRUE) 
par(mfrow=c(3,1))
par(mar=c(3.0, 5.5, 1.5, 2)) # bottom, left, top, right  

plot(h1, xlim=XLIMIT, ylim=c(0, 25), cex.axis=CEX.AXIS, main='', ylab='')
title(main = 'Accuracy', cex.main=CEX.MAIN, ylab='Frequency', cex.lab=CEX.LAB) 

plot(h2, xlim=XLIMIT, main = '', cex.axis=CEX.AXIS, ylab='') 
title(main = 'Average distance', cex.main=CEX.MAIN, ylab='Frequency', cex.lab=CEX.LAB)

plot(h3, xlim=XLIMIT, ylim=c(0, 55),main = '', cex.axis=CEX.AXIS, ylab='')
title(main = 'KL distance', cex.main=CEX.MAIN, ylab='Frequency', cex.lab=CEX.LAB)
dev.off()

# Histograms of circuit construction parameters 
# feature ratio cutoff, correlation cutoff, number of nodes
#---------------------------------------------------------- 
BREAKS <- 60
par(mfrow=c(3,1))
par(mar=c(2.0, 5.5, 1.5, 2)) # bottom, left, top, right 
h1 <- hist(circuit_metrics.sim$FeatureRatio, breaks = BREAKS, main = 'Feature ratio cutoff') 
h2 <- hist(circuit_metrics.sim$AbsCor, breaks = BREAKS, main = 'Correlation cutoff') 
h3 <- hist(circuit_metrics.sim$Nodes, breaks = BREAKS, main = 'Number of nodes') 
#hist(circuit_metrics.sim$Nodes, breaks = 20, main = 'Number of nodes') 

fname_fig <- paste0(outdir, 'histogram-parameters.pdf')
pdf(file=fname_fig, width=6, height=8, paper = "special", onefile = TRUE) 
par(mfrow=c(3,1))
par(mar=c(3.0, 5.5, 1.5, 2)) # bottom, left, top, right 

plot(h1, xlim=c(0.05, 0.2), ylim=c(0, 25) ,cex.axis=CEX.AXIS, main='', ylab='')
title(main = 'Feature ratio cutoff', cex.main=CEX.MAIN, ylab='Frequency', cex.lab=CEX.LAB) 

plot(h2, xlim = c(0, 1), ylim = c(0, 50) , main = '', cex.axis=CEX.AXIS, ylab='') 
title(main = 'Correlation cutoff', cex.main=CEX.MAIN, ylab='Frequency', cex.lab=CEX.LAB)

plot(h3, xlim = c(0, 60), ylim = c(0, 10), main = '', cex.axis=CEX.AXIS, ylab='')
title(main = 'Number of nodes', cex.main=CEX.MAIN, ylab='Frequency', cex.lab=CEX.LAB)
dev.off()

min(circuit_metrics.sim$Nodes)
sum(circuit_metrics.sim$Nodes<=5)


# Plot circuit construction parameters vs average distance
# feature ratio cutoff, correlation cutoff, number of nodes
#---------------------------------------------------------- 
fname_fig <- paste0(outdir, 'parameters-vs-average.dist.pdf')
pdf(file=fname_fig, width=6, height=8, paper = "special", onefile = TRUE) 
par(mfrow=c(3,1))
par(mar=c(5.0, 5.5, 1.5, 2)) # bottom, left, top, right 

plot(circuit_metrics.sim$FeatureRatio, circuit_metrics.sim$AvgDist, 
     main='', xlab='', ylab='',
     xlim=c(0.05, 0.2), ylim=c(0, 0.10), 
     cex.axis=CEX.AXIS)
title(xlab = 'Feature ratio cutoff', ylab='Average distance', cex.lab=CEX.LAB) 

plot(circuit_metrics.sim$Nodes, circuit_metrics.sim$AvgDist, 
     main='', xlab='', ylab='',
     xlim = c(0, 60), ylim = c(0, 0.10),
     cex.axis=CEX.AXIS)
title(xlab = 'Number of nodes', ylab='Average distance', cex.lab=CEX.LAB)

plot(circuit_metrics.sim$AbsCor, circuit_metrics.sim$AvgDist, 
     main='', xlab='', ylab='',
     xlim = c(0, 1.0), ylim = c(0, 0.10),
     cex.axis=CEX.AXIS)
title(xlab = 'Correlation cutoff', ylab='Average distance', cex.lab=CEX.LAB)
dev.off()


# Parameters vs metrics 
#======================
#Feature ratio
fname_fig <- paste0(outdir, 'featureRatio-vs-metrics.pdf')
pdf(file=fname_fig, width=6, height=8, paper = "special", onefile = TRUE) 
par(mfrow=c(3,1))
par(mar=c(5.0, 5.5, 1.5, 2)) # bottom, left, top, right 

plot(circuit_metrics.sim$FeatureRatio, circuit_metrics.sim$Accuracy, 
     main='', xlab='', ylab='',
     xlim=c(0.05, 0.2), #ylim=c(0, 0.10), 
     cex.axis=CEX.AXIS)
title(#xlab = 'Feature ratio cutoff', 
      ylab='Accuracy', cex.lab=CEX.LAB) 

plot(circuit_metrics.sim$FeatureRatio, circuit_metrics.sim$AvgDist, 
     main='', xlab='', ylab='',
     xlim=c(0.05, 0.2), #ylim = c(0, 0.10),
     cex.axis=CEX.AXIS)
title(#xlab = 'Feature ratio cutoff', 
      ylab='Average distance', cex.lab=CEX.LAB)

plot(circuit_metrics.sim$FeatureRatio, circuit_metrics.sim$flexibility, 
     main='', xlab='', ylab='',
     xlim=c(0.05, 0.2), #ylim = c(0, 0.10),
     cex.axis=CEX.AXIS)
title(xlab = 'Feature ratio cutoff', ylab='Flexibility', 
      cex.lab=CEX.LAB)
dev.off()


#Number of TFs
fname_fig <- paste0(outdir, 'TopTFs-vs-metrics.pdf')
pdf(file=fname_fig, width=6, height=8, paper = "special", onefile = TRUE) 
par(mfrow=c(3,1))
par(mar=c(5.0, 5.5, 1.5, 2)) # bottom, left, top, right 

plot(circuit_metrics.sim$TopTFs, circuit_metrics.sim$Accuracy, 
     main='', xlab='', ylab='',
     xlim=c(4, 32), 
     #ylim=c(0, 0.10), 
     cex.axis=CEX.AXIS)
title(#xlab = 'Feature ratio cutoff', 
  ylab='Accuracy', cex.lab=CEX.LAB) 

plot(circuit_metrics.sim$TopTFs, circuit_metrics.sim$AvgDist, 
     main='', xlab='', ylab='',
     xlim=c(4, 32), #ylim = c(0, 0.10),
     cex.axis=CEX.AXIS)
title(#xlab = 'Feature ratio cutoff', 
  ylab='Average distance', cex.lab=CEX.LAB)

plot(circuit_metrics.sim$TopTFs, circuit_metrics.sim$flexibility, 
     main='', xlab='', ylab='',
     xlim=c(4, 32), #ylim = c(0, 0.10),
     cex.axis=CEX.AXIS)
title(xlab = 'Top TFs', ylab='Flexibility', 
      cex.lab=CEX.LAB)
dev.off()


# Corr cut off
fname_fig <- paste0(outdir, 'AbsCorr-vs-metrics.pdf')
pdf(file=fname_fig, width=6, height=8, paper = "special", onefile = TRUE) 
par(mfrow=c(3,1))
par(mar=c(5.0, 5.5, 1.5, 2)) # bottom, left, top, right 

plot(circuit_metrics.sim$AbsCor, circuit_metrics.sim$Accuracy, 
     main='', xlab='', ylab='',
     xlim = c(0, 1.0), 
     #ylim=c(0, 0.10), 
     cex.axis=CEX.AXIS)
title(#xlab = 'Feature ratio cutoff', 
  ylab='Accuracy', cex.lab=CEX.LAB) 

plot(circuit_metrics.sim$AbsCor, circuit_metrics.sim$AvgDist, 
     main='', xlab='', ylab='',
     xlim = c(0, 1.0), #ylim = c(0, 0.10),
     cex.axis=CEX.AXIS)
title(#xlab = 'Feature ratio cutoff', 
  ylab='Average distance', cex.lab=CEX.LAB)

plot(circuit_metrics.sim$AbsCor, circuit_metrics.sim$flexibility, 
     main='', xlab='', ylab='',
     xlim = c(0, 1.0), #ylim = c(0, 0.10),
     cex.axis=CEX.AXIS)
title(xlab = 'Abs Correlation', ylab='Flexibility', 
      cex.lab=CEX.LAB)
dev.off()

unique(sort(circuit_metrics.sim$FeatureRatio))
unique(sort(circuit_metrics.sim$TopTFs))

# Select top circuits 
#====================
circuits.sele <- circuit_metrics.sim[circuit_metrics.sim$idxBoth<=13, ]  
dim(circuits.sele) # 5 circuits 

# select a subset of columns:
circuits.subset <- circuits.sele[, c("TopTFs", "FeatureRatio", "AbsCor", 
                                     "Nodes", "Interactions", "PosInt", 
                                     "Accuracy", "KLdist", "AvgDist")] 

dim(circuits.subset)
write.csv(circuits.subset, file = paste0(outdir, 'top.circuits.csv'), row.names = T, quote = F)

# load sRACIPE object for the selected circuits:
racipe.list <- list()
for(circuit_idx in rownames(circuits.subset)){
   print(circuit_idx) 
   racipe.list[[circuit_idx]] <- readRDS(file = paste('../circuits.sim/circuit_simulated_', 
                                  circuit_idx, '.rds', sep = '') ) 
   #break()
}

length(names(racipe.list))



# Find all the genes in all the selected circuits
#------------------------------------------------
nodes.list <- list()
for(circuit_idx in names(racipe.list)){ 
   racipe <- racipe.list[[circuit_idx]]
   nodes.list[[circuit_idx]] <- rownames(racipe) 
}

names(nodes.list)
nodes.all <- unlist(nodes.list)
length(unique(nodes.all)) # 51


# Explore one selected circuit
#------------------------------
names(racipe.list)
circuit_idx <- names(racipe.list)[1]  

# import simulated object:
racipe <- racipe.list[[circuit_idx]]

# plot circuit:
sRACIPE::sracipePlotCircuit(racipe, plotToFile = F)

# plot simulation data: 
sRACIPE::sracipePlotData(racipe, plotToFile = F) 

# 
fname.hS <- paste(hS.dir, 'hS_', circuit_idx, '.rds', sep = '')
fname.hS
hS <- readRDS(file = fname.hS)
hS$AvgDist
hS$simulated.cluster.freq[2] + hS$simulated.cluster.freq[3]
hS$KL

# plot reference data:
heatmap.2(hS$dataReference, trace = 'none') 

# plot simulated.refCor 
# hS$simulated.refCor: pearson correlation between data.REF and data.sim 
simulated.refCor <- hS$simulated.refCor
dim(simulated.refCor)

graphics::image(hS$simulated.refCor)
#graphics::image(t(hS$simulated.refCor))
