#!/usr/bin/env Rscript

remove(list = ls()) 

# library(sRACIPE)

NO_CLUSTERS.REF <- 2

XLIMIT <- c(-11, 5)
YLIMIT <- c(-4, 5)

# cluster_names <- c('CL0', 'CL1', 'CL2', 'CL3')
# names(cluster_names) <- c('Uncertain', 'AML', 'Hybrid', 'Untreated') 
# color_vector = c("black", "red", "dodgerblue", "orange") # CL0, CL1, CL2, CL3

cluster_names <- c('CL1', 'CL2')
names(cluster_names) <- c('Untreated', 'AML') 
color_vector = c("dodgerblue", "red") # CL1, CL2

figdir <- './figs.wt.clusters/'
dir.create(figdir)

# Load pca data 
#---------------
datadir <- './data/' 

NO_META_COLS = 1
fname_data <- '../racipe.wt/data/racipe.models.wt.labeled.csv'
mydata = read.csv(file=fname_data, header=TRUE)
dim(mydata)
mydata.saved <- mydata

# remove data in CL3 (NULL models)
mydata <- mydata[mydata$CLUSTER_NO!=3, ] 
dim(mydata)
ndata <- mydata[, (NO_META_COLS+1):ncol(mydata)] 
pca <- prcomp(ndata) 

# PC contributions 
#-----------------
eigs <- pca$sdev^2
stat <- rbind(SD = sqrt(eigs),
              Proportion = eigs/sum(eigs),
              Cumulative = cumsum(eigs)/sum(eigs))

# project the knockdown simulation data on pca from the wt simulations
#---------------------------------------------------------------------
prData.wt <- scale(mydata[, 2:ncol(mydata)] , pca$center, pca$scale) %*%  pca$rotation
mainstr <- paste(' Models: Untreated' , sep = '')

WIDTH <- 6  
HEIGHT <- 5 
figname <- paste(figdir, 'clusters_on_PCs.wt-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
plot(pca$x[,1], pca$x[,2], #xlim=XLIMIT, ylim=YLIMIT, 
     main=mainstr, col='gray', pch=19, cex=0.05, 
     xlab=paste('PC1 (', format(stat[2,1]*100, digits = 5), '%)' , sep = ''), 
     ylab=paste('PC2 (', format(stat[2,2]*100, digits = 4), '%)' , sep = '')
) 
for(cluster.no in 1:NO_CLUSTERS.REF){
  prCluster=prData.wt[mydata$CLUSTER_NO==cluster.no,]  
  if(is.vector(prCluster)) { 
    points(prCluster[1], prCluster[2], pch=19,cex=0.5, col=color_vector[cluster.no]) 
  }
  else{
    points(prCluster[,1], prCluster[,2],pch=19,cex=0.5, col=color_vector[cluster.no])  
  }
}
#legend(-1, (YLIMIT[1]+7), names(cluster_names), color_vector, bty='n') 
dev.off()

# Plot legend
#------------
WIDTH <- 3
HEIGHT <- 3
figname <- paste(figdir, 'legend-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special') 
plot(pca$x[,1], pca$x[,2], xlim=XLIMIT, ylim=YLIMIT, 
     col='white', pch=19, cex=0.005, xlab='PC1', ylab='PC2') 
legend(XLIMIT[1], (YLIMIT[2]-1), names(cluster_names), color_vector, bty='n') 
dev.off()
