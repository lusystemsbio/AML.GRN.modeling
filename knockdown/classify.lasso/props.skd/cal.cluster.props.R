remove(list = ls()) 

library(sRACIPE)
library(tidyverse)
library(glmnet)

NO_CLUSTERS.REF <- 2 # 3

outdir  <- './data/'
dir.create(outdir )

outdir.kd.clustered <- './data.sim.clustered/' 
dir.create(outdir.kd.clustered)

source('./functions.classify_models.R')

# Load trained model 
#-------------------
model <- readRDS('../racipe.wt/data/model.glmnet.rds')

# Load training data 
#-------------------
NO_META_COLS = 1
fname_data <- '../racipe.wt/data/racipe.models.wt.labeled.csv'
mydata = read.csv(file=fname_data, header=TRUE)
dim(mydata)

# Load mean and sd of models (with no knockdown) 
#-----------------------------------------------
fname_data <- '../racipe.wt/data/mean_sd.racipe.models.wt.csv'
mean_sd.wt.df = read.csv(file=fname_data, header=TRUE)

means.wt <- mean_sd.wt.df$mean
names(means.wt) <- mean_sd.wt.df$tf
sds.wt <- mean_sd.wt.df$sd
names(sds.wt) <- mean_sd.wt.df$tf
# sum(mean_sd.wt.df$tf==rn)

# Load kd simulations 
#--------------------
data.dir.sim <- '../../data.sim/skd/' 
fnames <- sort(list.files(data.dir.sim))
fnames
length(fnames)
fname <- fnames
xtract_gene_name <- function(fname){
   fname.suffix <- strsplit(fname, split = '-', fixed = TRUE)[[1]][2]
   gene_name <- strsplit(fname.suffix, split = '.', fixed = TRUE)[[1]][1]
   return(gene_name)
}

gene_names <- sapply(fnames, function(fname) xtract_gene_name(fname)) 
gene_names
length(gene_names)


#cluster_props <- as.data.frame(matrix(nrow = (length(fnames)+1), ncol = (NO_CLUSTERS.REF+1))) 
cluster_props <- as.data.frame(matrix(nrow = (length(fnames)+1), ncol = NO_CLUSTERS.REF)) 

rownames(cluster_props) <- c('Untreated', gene_names)
colnames(cluster_props) <- c('Cluster_1', 'Cluster_2')


fname <- fnames[1]
for(fname in fnames){
   print(fname)
   gene_name <- xtract_gene_name(fname)
   racipe <- readRDS(file = paste(data.dir.sim, fname, sep = ''))   
   
   geneExpression.kd <- assay(racipe,1) 
   geneExpression.kd <- log2(1+geneExpression.kd) 
   
   ndata.kd <- t(normalize_by_wt_mean_and_sd(geneExpression=geneExpression.kd, 
                                             means = means.wt, sds=sds.wt)) 
   
   probabilities <- model %>% predict(newx = ndata.kd, type='response') 
   predicted.classes <- ifelse(probabilities > 0.5, 1, 2) 
   
   # obtain cluster proportions: 
   cluster_props[gene_name,] <- c(sum(predicted.classes==1), 
                                  sum(predicted.classes==2)) 
   
   # save ndata.kd with cluster id attached
   ndata.kd.clustered <- cbind(predicted.classes, ndata.kd)  
   colnames(ndata.kd.clustered) <- c('CLUSTER_NO', colnames(ndata.kd))
   fname.ndata <- paste(outdir.kd.clustered, 'racipe-', gene_name, '.csv', sep = '')
   write.csv(ndata.kd.clustered, fname.ndata, quote = F, row.names = F)
   # break()
}

# Calculate cluster proportions for WT models using nnet preditor
#-----------------------------------------------------------------
probabilities <- model %>% predict(newx = as.matrix(mydata[, 2:ncol(mydata)]), type='response') 
predicted.classes <- ifelse(probabilities > 0.5, 1, 2) # WT models
cluster_props["Untreated",] <- c(sum(predicted.classes==1), 
                                 sum(predicted.classes==2))   

fname.out <- paste(outdir , 'cluster_props.csv', sep = '')
write.csv(cluster_props, file = fname.out, quote = F)
