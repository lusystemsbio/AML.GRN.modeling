#!/usr/bin/env Rscript

remove(list = ls())

sampleSizes <- c('01K', '02K', '03K', '04K', '05K', '06K', '07K', '08K', '09K', '10K')

meanSD.df <- as.data.frame(matrix(nrow = length(sampleSizes), ncol = 3))
colnames(meanSD.df) <- c('sample_size', 'mean', 'sd')  
rownames(meanSD.df) <- sampleSizes

for(sampleSize in sampleSizes){ 
   print(sampleSize)
   fname.input <- paste("./comb.metrics.SSIZE.", sampleSize, '/results/', 
                        'metric.comb.sorted.csv', sep = '') 
   print(fname.input)
   metric.comb <- read.csv(file = fname.input, row.names = 1) 
   
   meanSD.df[sampleSize, ] <- c(sampleSize, mean(metric.comb[,"sd"]), sd(metric.comb[,"sd"]))
   #break()
}


meanSD.df$sample_size <- factor(meanSD.df$sample_size, levels = sampleSizes)
meanSD.df$mean <- as.numeric(meanSD.df$mean)
meanSD.df$sd <- as.numeric(meanSD.df$sd)

library(ggplot2)

p1 <- ggplot(meanSD.df) + 
   geom_bar( aes(x=sample_size, y=mean), 
             stat="identity", fill="forestgreen", alpha=0.5) +
   geom_errorbar(aes(x=sample_size, 
                     ymin=mean-sd, ymax=mean+sd), width=0.4, 
                 colour="orange", alpha=0.9, size=1.5) +
   ggtitle("mean of SDs")

p1

WIDTH <- 5
HEIGHT <- 5
outdir <- './figs/'
dir.create(outdir)
fname.out <- paste0(outdir, 'mean.SDs.pdf')
ggsave(filename = fname.out, p1, width = WIDTH, height = HEIGHT)

 
#plot(meanSD.df$mean)

