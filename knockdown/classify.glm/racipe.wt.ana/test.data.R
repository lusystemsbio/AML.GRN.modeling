#!/usr/bin/env Rscript
remove(list = ls()) 

datadir <- '../racipe.wt/data/'

fname <- paste(datadir, 'racipe.models.wt.labeled.csv', sep='') 
mydata <- read.csv(fname)
mydata <- mydata[mydata$CLUSTER_NO!=3,]

fname.out <- paste('./data/sim.exp.data.csv', sep='')
write.csv(mydata, file = fname.out, quote = F, row.names = FALSE) 





hist(mydata$ATM)
hist(mydata$E2F1)
#plot(mydata$CLUSTER_NO, mydata$E2F1)
#plot(mydata$E2F1, mydata$CLUSTER_NO)


figdir <- './figs.clustered.data/' 
dir.create(figdir)
WIDTH <- 8
HEIGHT <- 12
figname <- paste(figdir, 'hist-genes-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
par(mfrow=c(5,6)) 
par(mar=c(2.5,1.1,1.1,1.1)) 
for(gene.name in colnames(mydata)[2:dim(mydata)[2]]){
  hist(mydata[,gene.name], xlab = '', ylab = '', main = gene.name)
   #break()
}
dev.off()


WIDTH <- 8
HEIGHT <- 12
figname <- paste(figdir, 'gene.exp-vs-cluster.no-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
par(mfrow=c(5,6)) 
par(mar=c(2.5,1.1,1.1,1.1))
for(gene.name in colnames(mydata)[2:dim(mydata)[2]]){
  plot(mydata[,gene.name], mydata$CLUSTER_NO, xlab = gene.name, ylab = 'clusters', 
       main = gene.name)
  #break()
}
dev.off()


WIDTH <- 8
HEIGHT <- 12
figname <- paste(figdir, 'cluster.no-vs-gene.exp-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
par(mfrow=c(5,6)) 
par(mar=c(2.5,1.1,1.1,1.1))
for(gene.name in colnames(mydata)[2:dim(mydata)[2]]){
  plot(mydata$CLUSTER_NO, mydata[,gene.name], xlab = gene.name, ylab = 'clusters', 
       main = gene.name)
  #break()
}
dev.off()
