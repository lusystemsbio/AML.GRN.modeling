#!/usr/bin/env Rscript

remove(list = ls()) 

datadir <- './data/'

fname <- paste(datadir, 'geneExpression.zscore.csv', sep='') 
mydata <- read.csv(fname)
#mydata <- mydata[mydata$CLUSTER_NO!=3,]

hist(mydata$ATM)
hist(mydata$E2F1)
#plot(mydata$CLUSTER_NO, mydata$E2F1)
#plot(mydata$E2F1, mydata$CLUSTER_NO)

figdir <- './figs.hist/' 
dir.create(figdir)
WIDTH <- 8
HEIGHT <- 12
figname <- paste(figdir, 'hist-genes-zscores-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
par(mfrow=c(5,6)) 
par(mar=c(2.5,1.1,1.1,1.1)) 
for(gene.name in colnames(mydata)){
  hist(mydata[,gene.name], xlab = '', ylab = '', main = gene.name)
  #break()
}
dev.off()
