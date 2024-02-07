
remove(list = ls()) 

library(sRACIPE)

datadir <- './data/'
dir.create(datadir)

circuit_idx <- '0.09-32-0.85'

# import simulated data
#----------------------
racipe <- readRDS(file = paste('../../../networks/circuits.sim/circuit_simulated_', 
                               circuit_idx, '.rds', sep = '') )

geneExpression <- assay(racipe,1) 
geneExpression <- log2(1+geneExpression)
means <- rowMeans(geneExpression)
sds <-  apply(geneExpression, 1, sd)

# Store mean and sd of the untreated data 
#-----------------------------------------
mean_sd.df <- cbind(names(means), means, sds)
colnames(mean_sd.df) <- c('tf' ,'mean', 'sd')

fname.out <- paste(datadir, 'mean_sd.racipe.models.wt.csv', sep='')
write.csv(mean_sd.df, file = fname.out, quote = F, row.names = FALSE)


# store racipe data (log transfer)
#----------------------------------
fname.out <- paste(datadir, 'geneExpression.csv', sep='')
write.csv(t(geneExpression), file = fname.out, quote = F, row.names = FALSE)

# standardize data
#---------------
mydata.df <- as.data.frame(t(geneExpression))

mean(mydata$ATM)
sd(mydata$ATM)
ndata.df <- mydata.df 
for(gene.name in rownames(mean_sd.df)){
  print(gene.name)
  ndata.df[, gene.name] <- (mydata.df[, gene.name] - as.numeric(mean_sd.df[gene.name, 'mean']))/as.numeric(mean_sd.df[gene.name, 'sd'])
}

fname.out <- paste(datadir, 'geneExpression.zscore.csv', sep='')
write.csv(ndata.df, file = fname.out, quote = F, row.names = FALSE)

# standardize data
#---------------
mydata.df <- as.data.frame(t(geneExpression))
ndata.df <- mydata.df 
for(gene.name in rownames(mean_sd.df)){
  print(gene.name)
  ndata.df[, gene.name] <- (mydata.df[, gene.name] - min(mydata.df[, gene.name]))/(max(mydata.df[, gene.name]) - min(mydata.df[, gene.name]))
}

fname.out <- paste(datadir, 'geneExpression.normed.csv', sep='')
write.csv(ndata.df, file = fname.out, quote = F, row.names = FALSE) 
