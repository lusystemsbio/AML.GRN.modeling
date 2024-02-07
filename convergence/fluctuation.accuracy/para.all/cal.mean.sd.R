#!/usr/bin/env Rscript
remove(list = ls()) 

# load accuracies - para all
fname <-'./results/accuracy.all.csv'
accuracy.all.df <- read.csv(file = fname, row.names = 1)
dim(accuracy.all.df)

mean.sd <- cbind(apply(accuracy.all.df, 1, mean), apply(accuracy.all.df, 1, sd))
colnames(mean.sd) <- c('mean', 'sd')

write.csv(format(mean.sd, digits = 3), file = "./results/accuracy.mean_sd.csv", row.names = T, quote = F) 

figdir <- './figs/'
dir.create(figdir)
WIDTH <- 8
HEIGHT <- 8
figname <- paste(figdir, 'accuracies-mean.sd-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')

par(mfcol = c(2, 1)) 
par(oma=c(3,1,3,3)) # b, l, t, r - all sides have 3 lines of space - outer margin
par(mar=c(1,4,1,1) + 0.1) # b, l, t, r - inner margin
plot(mean.sd[,'mean'], ylab = 'mean (accuracy)' )
plot(mean.sd[,'sd'], ylab='sd (accuracy)') 
dev.off()



