# Clean environment and set working direcotry
#--------------------------------------------
rm(list=ls())
setwd(getwd())

# Load library functions
#------------------------- 
libdir <- '../lib/'
#source(paste0(libdir, 'lib.aml.idh.R'))
source('./lib.data.tfs.R')

# Load TF list file
#------------------- 
fname_tfs <- './tfs.CTRL.IDH.csv'
tfs.df  <- read.csv(file = fname_tfs)


# Order TFS with Q val zero by Z.score
#=====================================
tfs.ordered  <- order.TFs_with_Q0(tfs.df)

# # # Save TFs after ordering
# #--------------------------- 
# fname_tfs.ordered  <- './tfs.CTRL.IDH.ordered.csv'
# write.csv(tfs.ordered, file = fname_tfs.ordered, row.names = FALSE, quote = FALSE)


#Save balanced TFs
#-----------------
fname.tfs.balanced <- './tfs.Netact.csv'
write.csv(format(tfs.ordered, digits=4), 
          file = fname.tfs.balanced, 
          quote = FALSE, row.names = F)


#Save selected total balanced tfs
#---------------------------------
fname.tfs <- './tfs.Netact.txt'
write.table(format(tfs.ordered, digits=4), 
            file = fname.tfs, quote = FALSE,
            row.names = FALSE, sep = '\t') 
