#!/usr/bin/env Rscript
remove(list = ls()) 

NO.TOP.CIRCUITS <- 10

circuit_metrics.sim <- read.csv(file = '../../rankings.mean_acc/results/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)

circuit_metrics.sim.list <- readRDS(file = paste0('../../networks/results/circuit_metrics.sim.list.rds'))
names(circuit_metrics.sim.list)

# load accuracies - para all
fname <-'../para.all/results/accuracy.all.csv'
accuracy.all.df <- read.csv(file = fname, row.names = 1)
dim(accuracy.all.df)


# create data structures
circuit_m.sim <- circuit_metrics.sim.list$S1

# check whetehr they are the same circuits
topCkts <- rownames(circuit_metrics.sim)[1:NO.TOP.CIRCUITS]
sum(rownames(circuit_m.sim) == topCkts)



figdir <- './figs/'
dir.create(figdir)
WIDTH <- 6  
HEIGHT <- 6  
figname <- paste(figdir, 'hist.accuracies-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special', onefile = TRUE)

par(mfcol = c(1, 1)) 
par(oma=c(3,3,3,3)) # b, l, t, r - all sides have 3 lines of space - outer margin
par(mar=c(1,1,4,1) + 0.1) # b, l, t, r - inner margin

idx <- rownames(accuracy.all.df)[1]

for(idx in rownames(circuit_metrics.sim)[1:NO.TOP.CIRCUITS]){
  # para all
  main.str <- paste(idx, ':', circuit_metrics.sim[idx, "Nodes"], sep = '')
  #hist(unlist(accuracy.all.df[idx,]), breaks = 50, main = main.str, xlim = xlimit) 
  hist(unlist(accuracy.all.df[idx,]), breaks = 50, main = main.str) 
  abline(v = mean(unlist(accuracy.all.df[idx,])), col="red", lwd=2, lty=1)
  #abline(v = circuit_metrics.sim[idx, "Accuracy"], col="red", lwd=2, lty=2) 
  
  
  for(sample_name in names(circuit_metrics.sim.list)){ 
    print(sample_name) 
    circuit_m.sim <- circuit_metrics.sim.list[[sample_name]]
    abline(v = circuit_m.sim[idx, "Accuracy"], col="blue", lwd=2, lty=2) 
  } 
  #break()
}

dev.off()





# flexibity 
#---------
fv <- c()
for(sn in names(circuit_metrics.sim.list)){
  fv <- c(fv, circuit_metrics.sim.list[[sn]]$flexibility)
}



figdir <- './figs/'
dir.create(figdir)
WIDTH <- 6  
HEIGHT <- 6  
figname <- paste(figdir, 'flexibility-bg-vs-sampled-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special', onefile = TRUE)

par(mfcol = c(1, 1)) 
par(oma=c(3,3,3,3)) # b, l, t, r - all sides have 3 lines of space - outer margin
par(mar=c(1,1,4,1) + 0.1) # b, l, t, r - inner margin

for(idx in rownames(circuit_metrics.sim)[1:NO.TOP.CIRCUITS]){
  # para all
  main.str <- paste(idx, ':', circuit_metrics.sim[idx, "Nodes"], sep = '')
  
  plot(fv, 1:length(fv), main = main.str)
  abline(v = circuit_metrics.sim[idx, "flexibility"], col="red", lwd=2, lty=1) 
  
  for(sample_name in names(circuit_metrics.sim.list)){ 
    print(sample_name) 
    circuit_m.sim <- circuit_metrics.sim.list[[sample_name]]
    abline(v = circuit_m.sim[idx, "flexibility"], col="blue", lwd=2, lty=2) 
  }
  #break()
}  
dev.off()



