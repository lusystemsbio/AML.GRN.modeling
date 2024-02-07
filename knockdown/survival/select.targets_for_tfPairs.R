#!/usr/bin/env Rscript
remove(list = ls()) 

circuit.idx <- '0.09-32-0.85'

datadir <- './data.target.tfPairs/'

# load tf pairs 
tfPairs <- readRDS(paste0(datadir, 'prodRate.pairs.rds'))  
class(tfPairs)
tfPairs$E2F4_TP53

# load TF-target DB
targetdb.list <- readRDS(file = '../../databases/targetDB.list.rds')
targetdb <- targetdb.list$`0.09`
names(targetdb)


tf.targets <- list()
tf.targetsthat_are_tfs <- list()
for(tfPair in names(tfPairs)){
  #print(tfPair)
  tf1 <- strsplit(tfPair, "_")[[1]][1]
  tf2 <- strsplit(tfPair, "_")[[1]][2] 
  targets.tfPair <- unique(c(tf1, tf2, targetdb[[tf1]], targetdb[[tf2]])) 
  print(length( targets.tfPair))
  tf.targets[[tfPair]] <- targets.tfPair 
  
  tf.targetsthat_are_tfs[[tfPair]] <- names(targetdb)[names(targetdb) %in% targets.tfPair] 
  #break()
}

names(tf.targets)
length(tf.targets$E2F4_TFDP1)
length(tf.targetsthat_are_tfs$E2F4_TFDP1)

length(names(tf.targets))
length(names(tf.targetsthat_are_tfs))


saveRDS(tf.targets, file = paste0(datadir, 'tf.targets.by_tfPair.rds'))
saveRDS(tf.targetsthat_are_tfs, file = paste0(datadir, 'tf.targets_that_are_tfs.bytfPair.rds'))

