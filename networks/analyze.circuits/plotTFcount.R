#!/usr/bin/env Rscript

remove(list = ls()) 

NO_TOP_CIRCUITS <- 10


tfCounts <- read.csv(file = '../results/tfCount.byMetod.csv', row.names = 1)

myColors <- c('black', 'red', 'blue', 'green')

plot(tfCounts$Nodes)
points(tfCounts$tfCount.Netact, col=myColors[2])
points(tfCounts$tfCount.MARINa, col=myColors[3])
points(tfCounts$tfCount.RI, col=myColors[4])


plot(tfCounts$tfCount.Netact, col=myColors[2], ylim = c(0, 30))
points(tfCounts$tfCount.MARINa, col=myColors[3])
points(tfCounts$tfCount.RI, col=myColors[4])




outdir <- './figs/'

