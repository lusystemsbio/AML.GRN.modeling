#!/usr/bin/env Rscript
remove(list = ls()) 

NO.TOP.CIRCUITS <- 10

# Load circuit summary
#-------------------- 
circuit_metrics.sim <- read.csv(file = '../rankingsBypercentile/results/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)
dim(circuit_metrics.sim)

# load accuracies
fname <-'./results.para/accuracy.dist.csv'
accuracy.df <- read.csv(file = fname, row.names = 1)
dim(accuracy.df)


perc.rank.3 <- function(x, xo)  (100-length(x[x <= xo])/length(x)*100)
#perc.rank.3 <- function(x, xo)  length(x[x >= xo])/length(x)*100

prank.samples.acc <- as.data.frame(matrix(nrow = nrow(accuracy.df), 
                                          ncol = ncol(accuracy.df)) ) 
rownames(prank.samples.acc) <- row.names(accuracy.df)
colnames(prank.samples.acc) <- colnames(accuracy.df)


prank.samples.flex <- as.data.frame(matrix(nrow = nrow(accuracy.df), 
                                           ncol = ncol(accuracy.df)))
rownames(prank.samples.flex) <- row.names(circuit_metrics.sim)[1:10]
colnames(prank.samples.flex) <- colnames(accuracy.df)

for(idx in rownames(accuracy.df)){ 
  print(idx) 
  prank.samples.acc[idx, ] <- sapply(unlist(accuracy.df[idx, ]), function (x) perc.rank.3(circuit_metrics.sim$Accuracy, x))
  #prank.samples.flex[sample_name] <- sapply(circuit_m.sim$flexibility, function (x) perc.rank.3(circuit_metrics.sim$flexibility, x))
  
}

figdir <- './figs/'
dir.create(figdir)
WIDTH <- 12 #8 #6 
HEIGHT <- 6 #12
figname <- paste(figdir, 'percentile.rank.accuracy-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')

par(mfrow = c(2, 5)) 
par(oma=c(3,3,3,3)) # b, l, t, r - all sides have 3 lines of space - outer margin
par(mar=c(1,1,4,1) + 0.1) # b, l, t, r - inner margin

xlimit <- c((min(prank.samples.acc)-0.01), c(max(prank.samples.acc)+0.01))

for(idx in rownames(prank.samples.acc)){
  print(idx)
  main.str <- paste(idx, ':', circuit_metrics.sim[idx, "Nodes"], sep = '')
  hist(unlist(prank.samples.acc[idx,]), breaks = 50, main = main.str) #, xlim = xlimit) 
  abline(v = mean(unlist(prank.samples.acc[idx,])), col="red", lwd=2, lty=1)
  abline(v = circuit_metrics.sim[idx, "idxAccuracy"], col="red", lwd=2, lty=2)
  #break()
}
dev.off()

