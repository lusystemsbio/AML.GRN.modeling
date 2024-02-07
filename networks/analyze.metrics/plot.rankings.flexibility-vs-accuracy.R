#!/usr/bin/env Rscript
remove(list = ls()) 

MINIMUM.CIRCUIT.SIZE <- 15

# Resources: 
# 1. http://r-statistics.co/Linear-Regression.html 

circuit_metrics.sim <- read.csv(file = '../results/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)
circuit_metrics.sim <- circuit_metrics.sim[circuit_metrics.sim$Nodes>=MINIMUM.CIRCUIT.SIZE, ]

figdir <- './figs.flexibility.accuracy/'
dir.create(figdir)


# accuracy vs flexibility  
#-------------------------
WIDTH <- 9
HEIGHT <- 9
fname_fig <- paste0(figdir, 'accuracy-vs-flexibility-', WIDTH, 'x', HEIGHT,'.pdf')
pdf(file=fname_fig, width=WIDTH, height=HEIGHT, paper = "special", onefile = TRUE) 
par(mfrow=c(1,1))
#par(mar=c(5.0, 5.5, 1.5, 2)) # bottom, left, top, right
plot(circuit_metrics.sim$idxAccuracy, circuit_metrics.sim$idxFlexibility, 
     xlab='accuracy ranking', ylab='flexibility ranking', 
     cex.axis = 2.0)
points(circuit_metrics.sim$idxAccuracy[1], circuit_metrics.sim$idxFlexibility[1], 
       cex=3.0, col='red', pch=19)
dev.off()



# flexibility vs accuracy
#-------------------------
WIDTH <- 9
HEIGHT <- 9
fname_fig <- paste0(figdir, 'flexibility-vs-accuracy-', WIDTH, 'x', HEIGHT,'.pdf')
pdf(file=fname_fig, width=WIDTH, height=HEIGHT, paper = "special", onefile = TRUE) 
par(mfrow=c(1,1))
#par(mar=c(5.0, 5.5, 1.5, 2)) # bottom, left, top, right
plot(circuit_metrics.sim$idxFlexibility, circuit_metrics.sim$idxAccuracy, 
     xlab='flexibility ranking', ylab='accuracy ranking')
points(circuit_metrics.sim$idxFlexibility[1], circuit_metrics.sim$idxAccuracy[1], 
       cex=3.0, col='red', pch=19) 
dev.off()




