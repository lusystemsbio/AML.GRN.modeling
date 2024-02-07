#!/usr/bin/env Rscript
remove(list = ls()) 

circuit.idx <- '0.09-32-0.85'

# top TFs from single and double KD analysis
#tfs <- c('RB1', 'STAT3', 'POU2F1', 'E2F1', 'MYC', 'TP53', 'E2F4', 'TFDP1')
tfs <- c('STAT3', 'RB1', 'POU2F1', 'MYCN', 'ETS2', 
         'NFIC', 'MYC',  'TP53', 'E2F4', 'TFDP1') 

length(tfs)  

# load TF-target DB
targetdb.list <- readRDS(file = '../../databases/targetDB.list.rds')
targetdb <- targetdb.list$`0.09`
names(targetdb)

tfs %in% names(targetdb)

targets <- c()
for(tf in tfs){
  targets <- c(targets, targetdb[[tf]])
}

targets
length(targets)
length(unique(targets))
tf.targets <- c(tfs, targets)
length(tf.targets)
length(unique(tf.targets))
tf.targets <- unique(tf.targets)
length(tf.targets) # 1876

tf.targets.that_are_tfs <- names(targetdb)[names(targetdb) %in% tf.targets]
length(tf.targets.that_are_tfs) # 179

length(tf.targets.that_are_tfs)
length(unique(tf.targets.that_are_tfs))
 
outdir <- './data.all.top.TFs/'
dir.create(outdir)
write.table(tf.targets, file = paste0(outdir, 'tf.targets.txt'), row.names = F, quote = F)
write.table(tf.targets.that_are_tfs, file = paste0(outdir, 'tf.targets.that_are_tfs.txt'), row.names = F, quote = F)
