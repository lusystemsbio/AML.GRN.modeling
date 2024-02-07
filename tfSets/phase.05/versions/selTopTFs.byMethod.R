#!/usr/bin/env Rscript
rm(list=ls()) 
setwd(getwd())

# Constants
#-----------
# number of TFs to be selected from each method:
#NO_TOP_TFS_FROM_EACH_METHOD.set <- seq(5, 50, 5)  
# set_1 <- seq(3, 15, 3)
# set_2 <- seq(15, 25, 2)
# set_2 <- set_2[2:length(set_2)]
NO_TOP_TFS_FROM_EACH_METHOD.set <- seq(4, 32, 4) 

# Output directories
#-------------------
outdir <- './data/'
dir.create(outdir)

# Import functions 
#-----------------
source('./functions.tfBalance.R')

# Load TFs (ordered) from each method
#--------------------------
fname.tfs.comb <- '../../data.tfs/TFs.combined.rds' 
TFs.combined <- readRDS(file = fname.tfs.comb)

# Load TF-target DB set
#---------------------------------------------------------------------------------------------- 
targetDB.list <- readRDS('../../databases/targetDB.list.rds') 

# Select TFs from each method and calculate the number of UP TFs in CTRL and TRTMT
#---------------------------------------------------------------------------------
coreTFs.list <- list()
for(fr in names(targetDB.list)){  
   targetDB.cur <- targetDB.list[[fr]] 
   for(top.TFs.count in NO_TOP_TFS_FROM_EACH_METHOD.set){ 
      # TFs from each method
      #tfs.cur <- as.character(TFs.combined[["Netact"]][,1])
      #tfsSele <- selectTOPtfs(TFs.combined=TFs.combined, targetDB=targetDB.cur, top.TFs.count=top.TFs.count)
      tfsSele <- selectTOPtfs.byMethod(TFs.combined=TFs.combined, 
                                       targetDB=targetDB.cur, 
                                       top.TFs.count=top.TFs.count) 
      coreTFs.list[[as.character(fr)]][[as.character(top.TFs.count)]] <- tfsSele  
      #break()
   }
   #break()
}

names(coreTFs.list)
length(names(coreTFs.list))

outdir <- './data/'
dir.create(outdir)
saveRDS(coreTFs.list, paste(outdir, 'coreTFs.rds', sep = ''))
