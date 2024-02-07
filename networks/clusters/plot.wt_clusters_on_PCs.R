remove(list = ls()) 

# best circuit
# XLIMIT <- c(-14, 6)
# YLIMIT <- c(-5, 5)

# largest circuit
XLIMIT <- c(-25, 6)
YLIMIT <- c(-6, 3)

NUM_CLUSTERS.REF <- 3
cluster_names <- c('CL1', 'CL2', 'CL3')
names(cluster_names) <- c('Untreated', 'AML', 'Hybrid')

# Define colors for each cluster 
#color_vector=c("orange", "red", "dodgerblue" , "magenta") 
#color_vector=c("deeppink", "darkgray", "dodgerblue") 
# color_vector=c("deeppink", "dodgerblue", "orange") 

# for best and largest circuits
color_vector=c("orange", "red",  "dodgerblue")
names(color_vector) <- names(cluster_names)

# for smallest circuit
names(cluster_names) <- c('AML', 'Untreated', 'Hybrid')
color_vector=c("red", "orange",   "dodgerblue")
names(color_vector) <- names(cluster_names)

library(sRACIPE)
source('./functions.classify_models.R')

outdir <- './data/'
dir.create(outdir)
figdir <- './figs.wt.clusters/'
dir.create(figdir)


# Load WT racipe models
#-----------------------
#circuit_idx <- "0.07-32-0.75"  # best network  
#circuit_idx <- "0.08-60-0.3"  #largest network 
circuit_idx <- "0.05-4-0.85"  #largest network 

fname <- paste('../circuits.sim/circuit_simulated_', circuit_idx, '.rds', sep = '')
#fname <- './circuit_simulated.rds'
racipe <- readRDS(file = fname)

geneExpression <- assay(racipe,1) 
geneExpression <- log2(1+geneExpression)
means <- rowMeans(geneExpression)
sds <-  apply(geneExpression, 1, sd)

ndata <- normalize_by_wt_mean_and_sd(geneExpression=geneExpression, 
                                     means=means, 
                                     sds=sds)
pca <- prcomp(t(ndata))

clusterCut <- obtain_clusterCut(ndata=ndata, num_clusters=NUM_CLUSTERS.REF)

# best network
# sum(clusterCut==1)  # 466 - Untreated 
# sum(clusterCut==2)  # 1415 - AML
# sum(clusterCut==3)  # 119 -  Hybrid 

# largest network
# sum(clusterCut==1)  # 220 - Untreated 
# sum(clusterCut==2)  # 1481 - AML
# sum(clusterCut==3)  # 299 -  Hybrid 

# smallest network
sum(clusterCut==1)  # 1225 - AML
sum(clusterCut==2)  # 687 - Untreated
sum(clusterCut==3)  # 88 -  Hybrid 


sum(clusterCut==1)*100/length(clusterCut)  
sum(clusterCut==2)*100/length(clusterCut)  
sum(clusterCut==3)*100/length(clusterCut) 

# Separate data into clusters 
#----------------------------
data_by_cluster <- obtain_data_by_cluster(data.sim=ndata, 
                                          clusterCut=clusterCut)
names(data_by_cluster)

# Obtain clusters projected on PCs
#---------------------------------
prData_by_cluster <- list() 
for(cluster_idx in names(data_by_cluster)){
   prData_by_cluster[[cluster_idx]] <- scale(t(data_by_cluster[[cluster_idx]]), pca$center, pca$scale) %*%  pca$rotation   #by using scale
}

# Save pca and projected wild type clusters 
#------------------------------------------
# pca.data <- list()
# pca.data[['pca']] <- pca 
# pca.data[['prData_by_cluster.wt']] <- prData_by_cluster
# fname.out <- paste(outdir, 'pca.data.rds', sep = '')
# saveRDS(pca.data, file = fname.out)

# Annotate cluster number with cluster names
#---------------------------------------------
names(prData_by_cluster) <- names(cluster_names)

WIDTH <- 7  
HEIGHT <- 5 
figname <- paste(figdir, 'clusters_from_data.sim.wt-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')

plot(pca$x[,1], pca$x[,2], xlim=XLIMIT, ylim=YLIMIT, 
     # main='data: models (Untreated)   dist: Euclidean',
     main='Models: untreated',
     col='gray', pch=19, cex=0.005, xlab='PC1', ylab='PC2') 

library(car) 
for(cluster_idx in names(prData_by_cluster)){ 
   prCluster=prData_by_cluster[[cluster_idx]]
   points(prCluster[,1], prCluster[,2],pch=19,cex=0.5, col=color_vector[cluster_idx])
   # dataEllipse(prData_by_cluster[[cluster_idx]][,1], prData_by_cluster[[cluster_idx]][,2],
   #             levels=c(0.8), add=TRUE,col=color_vector[cluster_idx],
   #             center.pch = 19,center.cex = 1.5, pch=1, cex=0.05, plot.points=TRUE,
   #             fill=TRUE, fill.alpha=0.1)
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
legend(XLIMIT[1], (YLIMIT[2]-3), names(cluster_names), color_vector, bty='n')
dev.off()

