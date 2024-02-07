#!/usr/bin/env Rscript
remove(list = ls()) 

NO.TOP.CIRCUITS <- 10

# Load circuit summary
#-------------------- 
datadir <- '../../rankingsBypercentile/results/'
circuit_metrics.sim <- read.csv(file = paste0(datadir , "./summary.circuits.sim.sortedByAcc_flex.csv"), row.names = 1)
circuit_metrics.sim <- circuit_metrics.sim[1:NO.TOP.CIRCUITS, ]

# load accuracies
fname <-'./results/accuracy.all.csv'
accuracy.df <- read.csv(file = fname, row.names = 1)
dim(accuracy.df)

idx <- rownames(accuracy.df)[1]
hist(unlist(accuracy.df[1,]), breaks = 50, main = idx)
abline(v = mean(unlist(accuracy.df[1,])), col="red", lwd=2, lty=1)
abline(v = circuit_metrics.sim[idx, "Accuracy"], col="red", lwd=2, lty=2)

# hist(unlist(accuracy.df[2,]), breaks = 50)
# hist(unlist(accuracy.df[3,]), breaks = 50)

figdir <- './figs/'
dir.create(figdir)
WIDTH <- 12 #8 #6 
HEIGHT <- 6 #12
figname <- paste(figdir, 'hist.accuracies-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
#par(mfrow = c(5, 2))   
par(mfrow = c(2, 5)) 
par(oma=c(3,3,3,3)) # b, l, t, r - all sides have 3 lines of space - outer margin
par(mar=c(1,1,4,1) + 0.1) # b, l, t, r - inner margin

xlimit <- c((min(accuracy.df)-0.01), c(max(accuracy.df)+0.01))
#xlimit <- c(0.90, 0.96)

for(idx in rownames(accuracy.df)){
  print(idx)
  main.str <- paste(idx, ':', circuit_metrics.sim[idx, "Nodes"], sep = '')
  hist(unlist(accuracy.df[idx,]), breaks = 50, main = main.str, xlim = xlimit) 
  abline(v = mean(unlist(accuracy.df[idx,])), col="red", lwd=2, lty=1)
  abline(v = circuit_metrics.sim[idx, "Accuracy"], col="red", lwd=2, lty=2)
  #break()
}
dev.off()


