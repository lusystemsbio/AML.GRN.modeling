
remove(list = ls()) 

library(sRACIPE)

NO_CLUSTERS.REF <- 2
#PROB_DIFF_THRESHOLD <- 0.05
# XLIMIT <- c(-6, 14)
# YLIMIT <- c(-10, 8)
XLIMIT <- c(-11, 5)
YLIMIT <- c(-4, 5)

# cluster_names <- c('CL0', 'CL1', 'CL2', 'CL3')
# names(cluster_names) <- c('Uncertain', 'AML', 'Hybrid', 'Unperturbed') 
# color_vector = c("black", "red", "dodgerblue", "orange") # CL0, CL1, CL2, CL3

cluster_names <- c('CL1', 'CL2')
names(cluster_names) <- c('Unperturbed', 'AML') 
color_vector = c("dodgerblue", "red") # CL1, CL2

figdir <- './figs.kd.clusters/'
dir.create(figdir)
source('./functions.classify_models.R')

# Load pca data 
#---------------
datadir <- './data/' 
#fname.in <- paste(datadir, 'pca.data.rds', sep = '')
#pca.data <- readRDS(file = fname.in)
#pca <- pca.data$pca
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

# Load clustered knockdown data 
#-----------------------------
datadir.sim <- './data.sim.clustered/' 
fnames <- sort(list.files(datadir.sim))
fnames
length(fnames)

fname <- fnames[1]
xtract_gene_name <- function(fname){
   fname.suffix <- strsplit(fname, split = '-', fixed = TRUE)[[1]][2]
   gene_name <- strsplit(fname.suffix, split = '.', fixed = TRUE)[[1]][1]
   return(gene_name)
}

gene_names <- sapply(fnames, function(fname) xtract_gene_name(fname)) 
gene_names
length(gene_names)

library(car) 

# Plot legend
#------------
WIDTH <- 12
HEIGHT <- 6
figname <- paste(figdir, 'legend-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')

plot(pca$x[,1], pca$x[,2], xlim=XLIMIT, ylim=YLIMIT, 
     col='white', pch=19, cex=0.005, xlab='PC1', ylab='PC2') 
#legend(XLIMIT[1], (YLIMIT[2]-1), names(cluster_names), color_vector, bty='n') 

groups <- c("Normal", "AML", "Background (Unperturbed)")
color_vector <- c(color_vector, "gray82")
legend(XLIMIT[1], (YLIMIT[2]-1), groups, color_vector, bty='n', horiz=TRUE) 
dev.off()
