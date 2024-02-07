#!/usr/bin/env Rscript

remove(list = ls()) 
circuit_metrics.sim <- read.csv(file = '../results/summary.circuits.sim.sorted.csv', row.names = 1)

max(circuit_metrics.sim$Nodes)

outdir <- './data/'
dir.create(outdir)

figdir <- './figs.bc/'
dir.create(figdir)

circuit_metrics.sim[1:10, 1:6]
unique(sort(circuit_metrics.sim[1:10, 4]))
length(unique(sort(circuit_metrics.sim[1:10, 4])))

circuit_metrics.sim[219:228, 1:6]
unique(sort(circuit_metrics.sim[219:228, 4]))
length(unique(sort(circuit_metrics.sim[219:228, 4])))


WIDTH <- 6 
HEIGHT <- 10 
figname <- paste(figdir, 'compare.indices.pdf', sep = '')
pdf(file = figname, paper = 'special', width = WIDTH, height = HEIGHT)
par(mfrow=c(3,1))
par(mar=c(1.1,4.0,2.1,1.1)) # bottom, left, top, right
plot(circuit_metrics.sim$Nodes, circuit_metrics.sim$idxAccuracy, #ylim=c(0, 1.0),
     xlab='', ylab='Accuracy')
plot(circuit_metrics.sim$Nodes, circuit_metrics.sim$idxAvgDist, #ylim=c(0, 1.0),
     xlab='', ylab='Avg Dist') 
par(mar=c(4.1,4.0,2.1,1.1)) # bottom, left, top, right
plot(circuit_metrics.sim$Nodes, circuit_metrics.sim$idxFlexibility, #ylim=c(0, 1.0),
     xlab='Nodes', ylab='Flexibility')
dev.off()



# Compare index (sorting index) pairs
#------------------------------------- 
WIDTH  <- 12
HEIGHT <- 8
figname <- paste(figdir, 'compare.index_pairs.pdf', sep = '')
pdf(file = figname, paper = 'special', width = WIDTH, height = HEIGHT)
par(mfrow=c(2,2))
par(mar=c(2.1,4.0,2.1,1.1)) # bottom, left, top, right
plot(circuit_metrics.sim$Nodes, circuit_metrics.sim$idxBoth,  
     xlab='', ylab='Accuracy + Avg Dist')
plot(circuit_metrics.sim$Nodes, circuit_metrics.sim$idxAccuracy+circuit_metrics.sim$idxFlexibility, 
     xlab='Nodes', ylab='Accuracy + Flexibility')
plot(circuit_metrics.sim$Nodes, circuit_metrics.sim$idxAvgDist+circuit_metrics.sim$idxFlexibility, 
     xlab='Nodes', ylab='Avg Dist + Flexibility')
#par(mar=c(4.1,4.0,2.1,1.1)) # bottom, left, top, right
plot(circuit_metrics.sim$Nodes, circuit_metrics.sim$idxTrio,  
     xlab='', ylab='Accuracy + Avg Dist + Flexibility') 
dev.off()




#=============================================== 
par(mfrow=c(1,1))
par(mar=c(4.1,4.0,2.1,1.1)) # bottom, left, top, right
plot(circuit_metrics.sim$Nodes, (circuit_metrics.sim$Accuracy)*1, 
     col='red', xlab='Nodes', ylab='metric', ylim=c(0, 1.0))
points(circuit_metrics.sim$Nodes, (circuit_metrics.sim$flexibility)*3, 
       col='black')

circuit_metrics.sim.tmp <- circuit_metrics.sim[circuit_metrics.sim$Nodes>=19 &
                                                  circuit_metrics.sim$Nodes<=22, ] 
plot(circuit_metrics.sim.tmp$Nodes, (circuit_metrics.sim.tmp$Accuracy)*1, 
     col='red', xlab='Nodes', ylab='metric', ylim=c(0, 1.0))
points(circuit_metrics.sim.tmp$Nodes, (circuit_metrics.sim.tmp$flexibility)*4, 
       col='black')


circuit_metrics.sim.19 <- circuit_metrics.sim.tmp[circuit_metrics.sim.tmp$Nodes==19, ] 
circuit_metrics.sim.20 <- circuit_metrics.sim.tmp[circuit_metrics.sim.tmp$Nodes==20, ]
circuit_metrics.sim.21 <- circuit_metrics.sim.tmp[circuit_metrics.sim.tmp$Nodes==21, ]
circuit_metrics.sim.22 <- circuit_metrics.sim.tmp[circuit_metrics.sim.tmp$Nodes==22, ]

outdir
fname.out <- paste(outdir, 'circuit_metrics.sim.19', '.csv', sep = '')
write.csv(circuit_metrics.sim.19, file = fname.out, quote = F, row.names = TRUE)

fname.out <- paste(outdir, 'circuit_metrics.sim.20', '.csv', sep = '')
write.csv(circuit_metrics.sim.20, file = fname.out, quote = F, row.names = TRUE)

fname.out <- paste(outdir, 'circuit_metrics.sim.21', '.csv', sep = '')
write.csv(circuit_metrics.sim.21, file = fname.out, quote = F, row.names = TRUE)

fname.out <- paste(outdir, 'circuit_metrics.sim.22', '.csv', sep = '')
write.csv(circuit_metrics.sim.22, file = fname.out, quote = F, row.names = TRUE)


#=============================================== 
# Transform into negative log10
# FlexAccuracy: flexibility*3 + Accuracy 
# circuits with low flexibility are NOT removed
#===============================================

# Remove circuits with nodes fewer than 10
# tmp.df <- circuit_metrics.sim[circuit_metrics.sim$Nodes <10, ]
# circuit_metrics.sim <- circuit_metrics.sim[circuit_metrics.sim$Nodes >=10, ]


WIDTH <- 6
HEIGHT <- 10
figname <- paste(figdir, 'log.metrics.pdf', sep = '')
pdf(file = figname, paper = 'special', width = WIDTH, height = HEIGHT)
par(mfrow=c(3,1))
par(mar=c(1.1,4.0,2.1,1.1)) # bottom, left, top, right

plot(circuit_metrics.sim$Nodes, -log10(circuit_metrics.sim$Accuracy), col='red', ylim=c(0, 2.0),
     xlab='', ylab='-log(metric)') 
points(circuit_metrics.sim$Nodes, -log10(circuit_metrics.sim$flexibility), col='black')

plot(circuit_metrics.sim$Nodes, -log10(circuit_metrics.sim$Accuracy)*2, col='red',
     xlab='Nodes', ylab='-log(metric)', ylim=c(0, 2.0))
points(circuit_metrics.sim$Nodes, -log10(circuit_metrics.sim$flexibility)/2, col='black')

par(mar=c(4.1,4.0,2.1,1.1)) # bottom, left, top, right
plot(circuit_metrics.sim$Nodes, -log10(circuit_metrics.sim$Accuracy)*2, col='red',
     xlab='Nodes', ylab='-log(metric)', ylim=c(0, 2.5))
points(circuit_metrics.sim$Nodes, -log10(circuit_metrics.sim$flexibility)/2, col='black')
points(circuit_metrics.sim$Nodes, -log10(circuit_metrics.sim$flexibility)/2
       -log10(circuit_metrics.sim$Accuracy)*2, col='blue')
dev.off()


# log transform
nlogFlex <- -log10(circuit_metrics.sim$flexibility)  
nlogAcc <- -log10(circuit_metrics.sim$Accuracy)

mean(nlogFlex) # 1.634605
median(nlogFlex) # 1.687805
sd(nlogFlex) # 0.2195559
mean(nlogAcc) # 0.099144
median(nlogAcc) # 0.05384261
sd(nlogAcc) # 0.1261049

# rescale
nlogFlex.z <- (nlogFlex-mean(nlogFlex))/sd(nlogFlex) 
nlogAcc.z <- (nlogAcc-mean(nlogAcc))/sd(nlogAcc) 
nlogComb.z <- nlogFlex.z + nlogAcc.z 

mydata <- cbind(rownames(circuit_metrics.sim), circuit_metrics.sim$Nodes, 
                nlogFlex.z, nlogAcc.z, nlogComb.z)
colnames(mydata) <- c('circuit_id', 'Nodes', 'nlogFlex', 'nlogAcc', 'nlogComb')
mydata <- as.data.frame(mydata) 
mydata$Nodes <- as.numeric(as.character(mydata$Nodes))
mydata$nlogAcc <- as.numeric(as.character(mydata$nlogAcc))
mydata$nlogFlex <- as.numeric(as.character(mydata$nlogFlex))
mydata$nlogComb <- as.numeric(as.character(mydata$nlogComb))

max(mydata$Nodes)


XLIMIT <- c((min(mydata$Nodes)-5), (max(mydata$Nodes)+5))

WIDTH <- 10
HEIGHT <- 6
figname <- paste(figdir, 'log.metrics.scaled.pdf', sep = '')
pdf(file = figname, paper = 'special', width = WIDTH, height = HEIGHT)
par(mfrow=c(1,1))
par(mar=c(4.1,4.0,2.1,1.1)) # bottom, left, top, right
plot(1, #mydata$Nodes, as.numeric(as.character(mydata$nlogAcc)), #col='red', 
     type="n", xlab='Nodes', xlim= XLIMIT, ylim=c(-6, 7)) 
points(mydata$Nodes, mydata$nlogAcc, col='red', pch = 18)
points(mydata$Nodes, mydata$nlogFlex, col='black', pch = 22)
points(mydata$Nodes, mydata$nlogComb, col='blue', pch=25)
legend(x=40, y=-1.5, 
       legend = c('negative log Accuracy', 'negative log Flexibility', 'Combined'), 
       col=c('red', 'black', 'blue'), pch=c(18, 22, 25))
dev.off()



WIDTH <- 4
HEIGHT <- 6
figname <- paste(figdir, 'hist.log.metrics.scaled.pdf', sep = '')
pdf(file = figname, paper = 'special', width = WIDTH, height = HEIGHT)
par(mfrow=c(3,1)) 
hist(mydata$nlogAcc, main = 'negative log Accuracy', xlab = '', xlim = c(-3, 4)) 
hist(mydata$nlogFlex, main = 'negative log Flexibility', xlab = '', xlim = c(-3, 4)) 
hist(mydata$nlogComb, main = 'Combined', xlab = '', xlim = c(-3, 4))
dev.off()

