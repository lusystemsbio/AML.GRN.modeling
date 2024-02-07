#!/usr/bin/env Rscript
remove(list = ls()) 

NO.TOP.CIRCUITS <- 10

# Load circuit summary
#-------------------- 
datadir <- '../../rankingsBypercentile/results/'
circuit_metrics.sim <- read.csv(file = paste0(datadir , "./summary.circuits.sim.sortedByAcc_flex.csv"), row.names = 1)
circuit_metrics.sim <- circuit_metrics.sim[1:NO.TOP.CIRCUITS, ]


# load accuracies - para 2
fname <-'../para.2/results/accuracy.all.csv'
accuracy.2.df <- read.csv(file = fname, row.names = 1)
dim(accuracy.2.df)

# load accuracies - para 3
fname <-'../para.3/results/accuracy.all.csv'
accuracy.3.df <- read.csv(file = fname, row.names = 1)
dim(accuracy.3.df)


# load accuracies - para 4
fname <-'./results/accuracy.all.csv'
accuracy.4.df <- read.csv(file = fname, row.names = 1)
dim(accuracy.4.df)


figdir <- './figs/'
dir.create(figdir)
WIDTH <- 12 #8 #6 
HEIGHT <- 6 #12
figname <- paste(figdir, 'hist.accuracies-', WIDTH, 'x', HEIGHT,'-top1-5-3runs.pdf', sep = '') 
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
#par(mfrow = c(5, 2))   
par(mfcol = c(3, 5)) 
par(oma=c(3,3,3,3)) # b, l, t, r - all sides have 3 lines of space - outer margin
par(mar=c(1,1,4,1) + 0.1) # b, l, t, r - inner margin

#xlimit <- c((min(accuracy.4.df)-0.01), c(max(accuracy.4.df)+0.01))
xlimit <- c(0.90, 0.96)

for(idx in rownames(accuracy.4.df)[1:5]){
  print(idx)
 
  # para 2
  main.str <- paste(idx, ':', circuit_metrics.sim[idx, "Nodes"], sep = '')
  hist(unlist(accuracy.2.df[idx,]), breaks = 50, main = main.str, xlim = xlimit) 
  abline(v = mean(unlist(accuracy.2.df[idx,])), col="red", lwd=2, lty=1)
  abline(v = circuit_metrics.sim[idx, "Accuracy"], col="red", lwd=2, lty=2)

  # para 3
  main.str <- paste(idx, ':', circuit_metrics.sim[idx, "Nodes"], sep = '')
  hist(unlist(accuracy.3.df[idx,]), breaks = 50, main = main.str, xlim = xlimit) 
  abline(v = mean(unlist(accuracy.3.df[idx,])), col="red", lwd=2, lty=1)
  abline(v = circuit_metrics.sim[idx, "Accuracy"], col="red", lwd=2, lty=2)
  
  # para 4
  main.str <- paste(idx, ':', circuit_metrics.sim[idx, "Nodes"], sep = '')
  hist(unlist(accuracy.4.df[idx,]), breaks = 50, main = main.str, xlim = xlimit) 
  abline(v = mean(unlist(accuracy.4.df[idx,])), col="red", lwd=2, lty=1)
  abline(v = circuit_metrics.sim[idx, "Accuracy"], col="red", lwd=2, lty=2)
  
  #break()
}
dev.off()


figdir <- './figs/'
dir.create(figdir)
WIDTH <- 12 #8 #6 
HEIGHT <- 6 #12
figname <- paste(figdir, 'hist.accuracies-', WIDTH, 'x', HEIGHT,'-top6-10-3runs.pdf', sep = '') 
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
#par(mfrow = c(5, 2))   
par(mfcol = c(3, 5)) 
par(oma=c(3,3,3,3)) # b, l, t, r - all sides have 3 lines of space - outer margin
par(mar=c(1,1,4,1) + 0.1) # b, l, t, r - inner margin

#xlimit <- c((min(accuracy.4.df)-0.01), c(max(accuracy.4.df)+0.01))
xlimit <- c(0.90, 0.96)

for(idx in rownames(accuracy.4.df)[6:10]){
  print(idx)
  
  # para 2
  main.str <- paste(idx, ':', circuit_metrics.sim[idx, "Nodes"], sep = '')
  hist(unlist(accuracy.2.df[idx,]), breaks = 50, main = main.str, xlim = xlimit) 
  abline(v = mean(unlist(accuracy.2.df[idx,])), col="red", lwd=2, lty=1)
  abline(v = circuit_metrics.sim[idx, "Accuracy"], col="red", lwd=2, lty=2)
  
  # para 3
  main.str <- paste(idx, ':', circuit_metrics.sim[idx, "Nodes"], sep = '')
  hist(unlist(accuracy.3.df[idx,]), breaks = 50, main = main.str, xlim = xlimit) 
  abline(v = mean(unlist(accuracy.3.df[idx,])), col="red", lwd=2, lty=1)
  abline(v = circuit_metrics.sim[idx, "Accuracy"], col="red", lwd=2, lty=2)
  
  # para 4
  main.str <- paste(idx, ':', circuit_metrics.sim[idx, "Nodes"], sep = '')
  hist(unlist(accuracy.4.df[idx,]), breaks = 50, main = main.str, xlim = xlimit) 
  abline(v = mean(unlist(accuracy.4.df[idx,])), col="red", lwd=2, lty=1)
  abline(v = circuit_metrics.sim[idx, "Accuracy"], col="red", lwd=2, lty=2)
  
  #break()
}
dev.off()
