#!/usr/bin/env Rscript

remove(list = ls()) 

datadir <- './data/'

fname <- paste(datadir, 'racipe.models.wt.labeled.csv', sep='') 
ndata.df <- read.csv(fname)
ndata.df <- ndata.df[ndata.df$CLUSTER_NO!=3,]

hist(ndata.df$ATM)
hist(ndata.df$E2F1)
plot(ndata.df$CLUSTER_NO, ndata.df$E2F1)
plot(ndata.df$E2F1, ndata.df$CLUSTER_NO)

figdir <- './figs.genes/' 
dir.create(figdir)
WIDTH <- 8
HEIGHT <- 12
figname <- paste(figdir, 'hist-genes-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
par(mfrow=c(5,6)) 
par(mar=c(2.5,1.1,1.1,1.1)) 
for(gene.name in colnames(ndata.df)[2:dim(ndata.df)[2]]){
  hist(ndata.df[,gene.name], xlab = '', ylab = '', main = gene.name)
  #break()
}
dev.off()


WIDTH <- 8
HEIGHT <- 12
figname <- paste(figdir, 'gene.exp-vs-cluster.no-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
par(mfrow=c(5,6)) 
par(mar=c(2.5,1.1,1.1,1.1))
for(gene.name in colnames(ndata.df)[2:dim(ndata.df)[2]]){
  plot(ndata.df[,gene.name], ndata.df$CLUSTER_NO, xlab = gene.name, ylab = 'clusters', 
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
for(gene.name in colnames(ndata.df)[2:dim(ndata.df)[2]]){
  plot(ndata.df$CLUSTER_NO, ndata.df[,gene.name], xlab = gene.name, ylab = 'clusters', 
       main = gene.name)
  #break()
}
dev.off()
