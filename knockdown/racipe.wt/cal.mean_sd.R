
remove(list = ls()) 

datadir <- './data/'
dir.create(datadir)

circuit_idx <- '0.09-32-0.85'

# import simulated data
#----------------------
racipe <- readRDS(file = paste('../../networks/circuits.sim/circuit_simulated_', 
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

