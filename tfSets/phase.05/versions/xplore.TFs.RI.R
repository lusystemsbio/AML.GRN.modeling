#!/usr/bin/env Rscript
rm(list=ls())
setwd(getwd())

# Load TFs from each method
#--------------------------
fname.tfs.comb <- '../../data.tfs/TFs.combined.rds' 
TFs.combined <- readRDS(file = fname.tfs.comb) 

# TFs deom RI method
#-------------------
tfs.RI <- TFs.combined$RI
colnames(tfs.RI)  

plot(tfs.RI$aggr.error )

# order: decreasing value of aggregate error:
tfs.RI.d <- tfs.RI[order(tfs.RI$aggr.error, decreasing = TRUE), ] 

# order: increasing value of aggregate error:
tfs.RI.i <- tfs.RI[order(tfs.RI$aggr.error, decreasing=FALSE), ]  


figdir <- './figs.RI/'
dir.create(figdir)
figname <- paste(figdir, 'RI-aggr.error.pdf', sep = '') 
WIDTH <- 5
HEIGHT <- 5 #10
pdf(file = figname, width=WIDTH, height=HEIGHT, paper='special')
#par(mfrow=c(2,1))
# plot(tfs.RI.d$aggr.error, col='gray', ylim = c(0, 0.015),
#      xlab = 'Regulators', ylab = 'Aggregate Error Change')
# abline(h=(mean(tfs.RI$aggr.error) + 1.5*sd(tfs.RI$aggr.error)))
plot(tfs.RI.i$aggr.error, col='gray', ylim = c(0, 0.015),
     xlab = 'Regulators', ylab = 'Aggregate Error Change')
abline(h=(mean(tfs.RI$aggr.error) + 1.5*sd(tfs.RI$aggr.error)))
dev.off()

fname.data <- paste(figdir, 'tfs.RI.csv', sep = '')
write.csv(tfs.RI.d, file = fname.data, quote = F, row.names = F)

TFs.xwen <- c('ERG', 'ELK1', 'RFX5', 'ELF1', 'PAX5', 'E2F1', 
              'DNMT1', 'MYC', 'E2F2', 'IRF9') 


TFs.xwen %in% as.character(tfs.RI.d$tf)
sum(as.character(tfs.RI.d$tf)=='DNMT1')

