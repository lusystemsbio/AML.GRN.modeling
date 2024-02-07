#!/usr/bin/env Rscript

remove(list = ls()) 

datadir <- './data/'

fname <- paste(datadir, 'mean_sd.racipe.models.wt.csv', sep='')
mean_sd.df <- read.csv(fname, row.names = 1)


fname <- paste(datadir, 'racipe.models.wt.labeled.csv', sep='') 
ndata.df <- read.csv(fname)

data.df <- as.data.frame(matrix(nrow = nrow(ndata.df), ncol = ncol(ndata.df)))
colnames(data.df) <- colnames(ndata.df)

col.no <- 1
for(gname in colnames(ndata.df)[2:ncol(ndata.df)]){
  print(gname)
  data.df[, gname] <- ndata.df[, gname] * mean_sd.df[gname, "sd"] + mean_sd.df[gname, "mean"] 
  #break()
}

data.df$CLUSTER_NO <- ndata.df$CLUSTER_NO


fname.out <- paste(datadir, 'racipe.models.unnormalized.csv', sep = '')
write.csv(data.df, fname.out, row.names = FALSE, quote = F)




