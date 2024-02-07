#!/usr/bin/env Rscript

remove(list = ls())
circuit_metrics.sim <- read.csv(file = '../results/summary.circuits.sim.sorted.csv', row.names = 1)
circuit_metrics.sim.s <- circuit_metrics.sim[order(circuit_metrics.sim$Accuracy, decreasing = TRUE), ]
write.csv(circuit_metrics.sim.s, './data/circuit_metrics.sortedByAcc.csv', 
          row.names = T, quote = F)

circuit_metrics.sim.s <- circuit_metrics.sim[order(circuit_metrics.sim$flexibility, decreasing = TRUE), ]
write.csv(circuit_metrics.sim.s, './data/circuit_metrics.sortedByFlex.csv', 
          row.names = T, quote = F)


flexibility.scaled <- (circuit_metrics.sim$flexibility-min(circuit_metrics.sim$flexibility))/(max(circuit_metrics.sim$flexibility)-min(circuit_metrics.sim$flexibility))
mean(circuit_metrics.sim$Accuracy)
median(circuit_metrics.sim$Accuracy)

mean(circuit_metrics.sim$flexibility)
median(circuit_metrics.sim$flexibility)

sd(circuit_metrics.sim$Accuracy)
sd(circuit_metrics.sim$flexibility)

min(circuit_metrics.sim.s$flexibility)
max(circuit_metrics.sim.s$flexibility)

mean_ratio <- mean(circuit_metrics.sim$Accuracy)/mean(circuit_metrics.sim$flexibility)
median_ratio <- median(circuit_metrics.sim$Accuracy)/median(circuit_metrics.sim$flexibility)

mean_ratio
median_ratio 


# Histograms of Accuracy and flexibility  
BREAKS <- 30
par(mfrow=c(2,2))
hist(circuit_metrics.sim$Accuracy, breaks = BREAKS, xlim = c(0, 1), 
     xlab = '', main = 'Accuracy') 
hist(circuit_metrics.sim$flexibility, breaks = BREAKS, xlim = c(0, 0.3), 
     xlab = '', ylab = '', main = 'Flexibility') 
hist(circuit_metrics.sim$flexibility*3, breaks = BREAKS, xlim = c(0, 1), 
     xlab = '', main = 'Flexibility') 
hist(flexibility.scaled, breaks = BREAKS, xlim = c(0, 1), 
     xlab = '', ylab = '', main = 'Flexibility') 


par(mfrow=c(2,2))
par(mar=c(2.1, 4.1, 1.1, 1.1)) # bottom, left, top, right
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$Accuracy, ylim=c(0, 1.0),
     xlab='', ylab='Accuracy')
plot(circuit_metrics.sim.s$Interactions, circuit_metrics.sim.s$Accuracy, ylim=c(0, 1.0),
     xlab='', ylab='')

par(mar=c(4.1, 4.1, 1.1, 1.1)) # bottom, left, top, right
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexibility, 
     ylim=c(0, 0.35), xlab='Nodes', ylab='flexibility')
# par(mar=c(2.1, 4.1, 1.1, 1.1)) # bottom, left, top, right
plot(circuit_metrics.sim.s$Interactions, circuit_metrics.sim.s$flexibility, 
     ylim=c(0, 0.35), xlab='Interactions', ylab='')





par(mfrow=c(2,2))
plot(circuit_metrics.sim$Accuracy, circuit_metrics.sim$flexibility, 
     xlab='Accuracy', ylab='Flexibility') 
plot(circuit_metrics.sim$Accuracy, circuit_metrics.sim$flexibility*3, 
     xlab='Accuracy', ylab='Flexibility x 3') 

tmp.df <- circuit_metrics.sim[circuit_metrics.sim$flexibility>=0.05, ]
dim(tmp.df)
plot(tmp.df$Accuracy, tmp.df$flexibility, xlab='Accuracy', ylab='Flexibility') 
plot(tmp.df$Accuracy, tmp.df$flexibility*3, xlab='Accuracy', ylab='Flexibility x 3') 


# remove circuits with fewer than 10 nodes 
tmp2.df <- tmp.df[tmp.df$Nodes>=15, ]
par(mfrow=c(2,1))
plot(tmp2.df$Accuracy, tmp2.df$flexibility, xlab='Accuracy', ylab='Flexibility') 
#plot(tmp2.df$Accuracy, tmp2.df$flexibility*3, xlab='Accuracy', ylab='Flexibility x 3')
plot(tmp2.df$Nodes, tmp2.df$flexibility, xlab='Nodes', ylab='Flexibility') 
#plot(tmp2.df$Nodes, tmp2.df$flexibility*3, xlab='Nodes', ylab='Flexibility x 3')


#========================================
# FlexAccuracy: flexibility  + Accuracy 
# circuits with low flexibility are NOT removed
#========================================
flexAccuracy <- circuit_metrics.sim$flexibility+circuit_metrics.sim$Accuracy

circuit_metrics.sim.2 <- cbind(circuit_metrics.sim, flexAccuracy)
# remove low flexibiliyt circuits
#circuit_metrics.sim.2 <- circuit_metrics.sim.2[circuit_metrics.sim.2$flexibility>=0.05, ]

circuit_metrics.sim.s <- circuit_metrics.sim.2[order(circuit_metrics.sim.2$flexAccuracy, 
                                                     decreasing = TRUE), ] 


plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$Accuracy, #ylim=c(0, 1.0),
     xlab='', ylab='Accuracy')
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexibility, #ylim=c(0, 1.0),
     xlab='', ylab='Flexibility')
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexibility*3, #ylim=c(0, 1.0),
     xlab='', ylab='Flexibility x 3')
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexAccuracy, #ylim=c(0, 1.0),
     xlab='Nodes', ylab='FlexAccuracy')




# Histograms of Accuracy and flexibility  
BREAKS <- 30
par(mfrow=c(2,2))
hist(circuit_metrics.sim.2$Accuracy, breaks = BREAKS, xlim = c(0, 1), 
     xlab = '', main = 'Accuracy') 
hist(circuit_metrics.sim.2$flexibility, breaks = BREAKS, xlim = c(0, 0.3), 
     xlab = '', ylab = '', main = 'Flexibility') 
hist(circuit_metrics.sim.2$flexibility*3, breaks = BREAKS, xlim = c(0, 1), 
     xlab = '', main = 'Flexibility x 3') 
hist(circuit_metrics.sim.2$flexAccuracy, breaks = BREAKS, xlim = c(0.0, 1.9), 
     xlab = '', ylab = '', main = 'FlexAccuracy')


par(mfrow=c(2,1))
par(mar=c(4.5,1.1,1.1,1.1)) # bottom, left, top, right
plot(-log10(circuit_metrics.sim$flexibility), ylim=c(0.0, 2.0), xlab='', ylab='')
points(-log10(circuit_metrics.sim$Accuracy), col='red') 

plot(-log10(circuit_metrics.sim$flexibility)/2, ylim=c(0.0, 2.0), xlab='Circuits', ylab='')
points(-log10(circuit_metrics.sim$Accuracy)*2, col='red') 


par(mfrow=c(2,1))
plot(circuit_metrics.sim$Nodes, -log10(circuit_metrics.sim$Accuracy), col='red', ylim=c(0, 2.0), 
     xlab='', ylab='-log(metric)')  
points(circuit_metrics.sim$Nodes, -log10(circuit_metrics.sim$flexibility), col='black')

plot(circuit_metrics.sim$Nodes, -log10(circuit_metrics.sim$Accuracy)*2, col='red', 
     xlab='Nodes', ylab='-log(metric)', ylim=c(0, 2.0))  
points(circuit_metrics.sim$Nodes, -log10(circuit_metrics.sim$flexibility)/2, col='black')



#========================================
# FlexAccuracy: flexibility + Accuracy 
# circuits with low flexibility are removed
#========================================
flexAccuracy <- circuit_metrics.sim$flexibility+circuit_metrics.sim$Accuracy

circuit_metrics.sim.2 <- cbind(circuit_metrics.sim, flexAccuracy)
# remove low flexibiliyt circuits
circuit_metrics.sim.2 <- circuit_metrics.sim.2[circuit_metrics.sim.2$flexibility>=0.05, ]

circuit_metrics.sim.s <- circuit_metrics.sim.2[order(circuit_metrics.sim.2$flexAccuracy, 
                                                     decreasing = TRUE), ] 

par(mfrow=c(2,2))
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$Accuracy, #ylim=c(0, 1.0),
     xlab='', ylab='Accuracy')
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexibility, #ylim=c(0, 1.0),
     xlab='', ylab='Flexibility')
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexibility*3, #ylim=c(0, 1.0),
     xlab='', ylab='Flexibility x 3')
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexAccuracy, #ylim=c(0, 1.0),
     xlab='Nodes', ylab='FlexAccuracy')


# Histograms of Accuracy and flexibility  
BREAKS <- 30
par(mfrow=c(2,2))
hist(circuit_metrics.sim.2$Accuracy, breaks = BREAKS, xlim = c(0, 1), 
     xlab = '', main = 'Accuracy') 
hist(circuit_metrics.sim.2$flexibility, breaks = BREAKS, xlim = c(0, 0.3), 
     xlab = '', ylab = '', main = 'Flexibility') 
hist(circuit_metrics.sim.2$flexibility*3, breaks = BREAKS, xlim = c(0, 1), 
     xlab = '', main = 'Flexibility x 3') 
hist(circuit_metrics.sim.2$flexAccuracy, breaks = BREAKS, xlim = c(0.0, 1.9), 
     xlab = '', ylab = '', main = 'FlexAccuracy')


#========================================
# FlexAccuracy: flexibility *3 + Accuracy 
# circuits with low flexibility are NOT removed
#========================================
flexAccuracy <- circuit_metrics.sim$flexibility*3+circuit_metrics.sim$Accuracy

circuit_metrics.sim.2 <- cbind(circuit_metrics.sim, flexAccuracy)
# remove low flexibiliyt circuits
#circuit_metrics.sim.2 <- circuit_metrics.sim.2[circuit_metrics.sim.2$flexibility>=0.05, ]

circuit_metrics.sim.s <- circuit_metrics.sim.2[order(circuit_metrics.sim.2$flexAccuracy, 
                                                     decreasing = TRUE), ] 

par(mfrow=c(2,2))
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$Accuracy, #ylim=c(0, 1.0),
     xlab='', ylab='Accuracy')
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexibility, #ylim=c(0, 1.0),
     xlab='', ylab='Flexibility')
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexibility*3, #ylim=c(0, 1.0),
     xlab='', ylab='Flexibility x 3')
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexAccuracy, #ylim=c(0, 1.0),
     xlab='Nodes', ylab='FlexAccuracy')


# Histograms of Accuracy and flexibility  
BREAKS <- 30
par(mfrow=c(2,2))
hist(circuit_metrics.sim.2$Accuracy, breaks = BREAKS, xlim = c(0, 1), 
     xlab = '', main = 'Accuracy') 
hist(circuit_metrics.sim.2$flexibility, breaks = BREAKS, xlim = c(0, 0.3), 
     xlab = '', ylab = '', main = 'Flexibility') 
hist(circuit_metrics.sim.2$flexibility*3, breaks = BREAKS, xlim = c(0, 1), 
     xlab = '', main = 'Flexibility x 3') 
hist(circuit_metrics.sim.2$flexAccuracy, breaks = BREAKS, xlim = c(0.0, 1.9), 
     xlab = '', ylab = '', main = 'FlexAccuracy')


#========================================
# FlexAccuracy: flexibility *3 + Accuracy 
# circuits with low flexibility are removed
#========================================
flexAccuracy <- circuit_metrics.sim$flexibility*3+circuit_metrics.sim$Accuracy

circuit_metrics.sim.2 <- cbind(circuit_metrics.sim, flexAccuracy)
# remove low flexibility circuits
circuit_metrics.sim.2 <- circuit_metrics.sim.2[circuit_metrics.sim.2$flexibility>=0.05, ]

circuit_metrics.sim.s <- circuit_metrics.sim.2[order(circuit_metrics.sim.2$flexAccuracy, 
                                                     decreasing = TRUE), ] 

par(mfrow=c(2,2))
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$Accuracy, #ylim=c(0, 1.0),
     xlab='', ylab='Accuracy')
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexibility, #ylim=c(0, 1.0),
     xlab='', ylab='Flexibility')
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexibility*3, #ylim=c(0, 1.0),
     xlab='', ylab='Flexibility x 3')
plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexAccuracy, #ylim=c(0, 1.0),
     xlab='Nodes', ylab='FlexAccuracy')
abline(h=0.42)


# Histograms of Accuracy and flexibility  
BREAKS <- 30
par(mfrow=c(2,2))
hist(circuit_metrics.sim.2$Accuracy, breaks = BREAKS, xlim = c(0, 1), 
     xlab = '', main = 'Accuracy') 
hist(circuit_metrics.sim.2$flexibility, breaks = BREAKS, xlim = c(0, 0.3), 
     xlab = '', ylab = '', main = 'Flexibility') 
hist(circuit_metrics.sim.2$flexibility*3, breaks = BREAKS, xlim = c(0, 1), 
     xlab = '', main = 'Flexibility x 3') 
hist(circuit_metrics.sim.2$flexAccuracy, breaks = BREAKS, xlim = c(0.3, 1.9), 
     xlab = '', ylab = '', main = 'FlexAccuracy')


plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$flexAccuracy, #ylim=c(0, 1.0),
     xlab='Nodes', ylab='FlexAccuracy')
abline(h=0.42)
circuits.sele <- circuit_metrics.sim.s[circuit_metrics.sim.s$Nodes>15 & 
                                          circuit_metrics.sim.s$Nodes < 20 & 
                                          circuit_metrics.sim.s$flexAccuracy <=0.42,  ]  

circuits.sele.2 <- circuits.sele[, c(1:6, 12:15, 18:19, 26)]

#========================
flexibility.scaled.s <- (circuit_metrics.sim.s$flexibility-min(circuit_metrics.sim.s$flexibility))/(max(circuit_metrics.sim.s$flexibility)-min(circuit_metrics.sim.s$flexibility))
par(mfrow=c(2,2))
plot(circuit_metrics.sim.s$Accuracy)
points(circuit_metrics.sim.s$flexibility*3, col='orange')
 
# plot(circuit_metrics.sim.s$Accuracy)
# points(flexibility.scaled, col='orange')

plot(circuit_metrics.sim.s$Nodes, circuit_metrics.sim.s$Interactions, 
     xlab='Nodes', ylab='Interactions')

plot(circuit_metrics.sim.s$Interactions/circuit_metrics.sim.s$Nodes)
plot(circuit_metrics.sim.s$Nodes)
plot(circuit_metrics.sim.s$Interactions)

subset.df <- circuit_metrics.sim.s[1:50, ]
plot(subset.df$Accuracy)
points(subset.df$flexibility*3, col='orange')

par(mfrow=c(3,1))

plot(subset.df$Nodes, xlab='')
plot(subset.df$Interactions, xlab='')
plot(subset.df$Interactions/subset.df$Nodes)

subset.df2 <- subset.df[10:20, ]

subset.df3 <- subset.df2[subset.df2$Accuracy>=0.6, ]

