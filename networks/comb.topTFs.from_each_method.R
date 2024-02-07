#!/usr/bin/env Rscript
rm(list=ls())
setwd(getwd())

NO_TOP_TFS_FROM_EACH_METHOD.set <- seq(5, 50, 5) 

source('./functions.R')

# Load TFs from each method
#--------------------------
fname.tfs.comb <- '../data.tfs/TFs.combined.rds' 
TFs.combined <- readRDS(file = fname.tfs.comb) 

# TFs from netact
#----------------
tfs.netact <- TFs.combined$NetAct
length(tfs.netact$tf)  

# order the data frame by increasing value of z and q values:  
tfs.netact <- tfs.netact[order(tfs.netact$zq, decreasing = FALSE), ]

# TFs from MARINa
#----------------
tfs.MARINa <- TFs.combined$MARINa
# order the data frame by increasing value of FDR 
tfs.MARINa  <- tfs.MARINa[order(tfs.MARINa$FDR, decreasing = FALSE), ] 


# TFs deom RI method
#-------------------
tfs.RI <- TFs.combined$RI
colnames(tfs.RI)  
# order the data frame by decreasing value of aggregate error:
tfs.RI <- tfs.RI[order(tfs.RI$aggr.error, decreasing = TRUE), ]

# Load combined TFs from three methods: NetAct, MARINa, and RI
#-------------------------------------------------------------
coreTFs.list <- list()
for(top.TFs.count in NO_TOP_TFS_FROM_EACH_METHOD.set){
   # obtain core TFs: 
   coreTFs.list[[as.character(top.TFs.count)]] <- sel.topTFs.from_each_method(TFs.combined=TFs.combined, 
                                                                              NO_TOP_TFS_FROM_EACH_METHOD=top.TFs.count) 
   # print(length(coreTFs))  
   # coreTFs.list[[as.character(top.TFs.count)]] <- coreTFs
}

names(coreTFs.list)
length(names(coreTFs.list))

saveRDS(coreTFs.list, './coreTFs.rds')

