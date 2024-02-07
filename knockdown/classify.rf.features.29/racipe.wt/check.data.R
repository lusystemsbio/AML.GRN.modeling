#!/usr/bin/env Rscript

remove(list = ls()) 

# check data - hS returned models vs simulated models 

circuit_idx <- '0.09-32-0.85'
NO_MODELS <- 10000


# import simulated data
#----------------------
racipe <- readRDS(file = paste('../../networks/circuits.sim/circuit_simulated_', 
                               circuit_idx, '.rds', sep = '') )

geneExpression <- assay(racipe,1) 
geneExpression <- log2(1+geneExpression)
means <- rowMeans(geneExpression)
sds <-  apply(geneExpression, 1, sd)

ndata <- normalize_by_wt_mean_and_sd(geneExpression=geneExpression, 
                                     means=means, 
                                     sds=sds)
pca.sim <- prcomp(t(ndata)) 

# import hS data
#---------------
hS.dir <- '../../networks/circuits.hS/' 
fname.hS <- paste(hS.dir, 'hS_', circuit_idx, '.rds', sep = '')
fname.hS
hS <- readRDS(file = fname.hS)
hS$AvgDist
hS$simulated.cluster.freq
hS$simulated.cluster.freq[2] + hS$simulated.cluster.freq[3]
hS$KL

data.sim.hS <- t(hS$dataSimulation)
dim(data.sim.hS)

pca.hS <- prcomp(data.sim.hS)


# project data
#-------------
plot(pca.hS$x[,1], pca.hS$x[,2], #xlim=XLIMIT, ylim=YLIMIT, 
     # main='data: models (Untreated)   dist: Euclidean',
     main='Models: untreated',
     col='gray', pch=19, cex=0.1, xlab='PC1', ylab='PC2')


plot(pca.sim$x[,1], pca.sim$x[,2], #xlim=XLIMIT, ylim=YLIMIT, 
     # main='data: models (Untreated)   dist: Euclidean',
     main='Models: untreated',
     col='gray', pch=19, cex=0.1, xlab='PC1', ylab='PC2')


cor(pca.hS$x[,1], pca.sim$x[,1])
cor(pca.hS$x[,2], pca.sim$x[,2]) 

cor(sort(pca.hS$x[,1]), sort(pca.sim$x[,1]))

plot(pca.hS$x[,1])
plot(pca.sim$x[,1]) 

plot(sort(pca.hS$x[,1]))
plot(sort(pca.sim$x[,1]) )




