# Clean environment and set working direcotry
#--------------------------------------------
rm(list=ls()) 
setwd(getwd())

# set directory paths
#-------------------
datadir_raw <- '../data.raw/'

#---------------------#
# CONSTANTS
#---------------------#
NO_PERMUTATIONS <- 10000
CONTROL <- 'CTRL'

#---------------------#
# Load library functions
#---------------------#
libdir <- '../lib/'
source('./lib.data.tfs.R')

#---------------------#
# Load expression data
#---------------------#
fname_data <- paste0(datadir_raw, "raw_brainarray.sele.txt") 
edata <- read.table(file = fname_data, 
                    header=TRUE, sep ='\t', 
                    row.names=1)
dim(edata)

#------------------------------------------------------#
# Create an ExpressionSet instance from expression data
#------------------------------------------------------#
eset.brain_array <- create.ExpressionSet_instance(edata)   

# Save eExpressionSet instance
#-----------------------------
save(eset.brain_array, file = paste('eset.brain_array.rda', sep = '', 
                                    collapse = NULL))

#--------------------#
# Calculate DE genes
#--------------------#
# create compare list
#compList = c("CTRL-IDH1", "CTRL-IDH2")  
# Calculate multi comparison DE genes:
# library(NetAct)
# de.results = MicroProcess(eset = eset.brain_array, compList = compList) 
# summary(de.results)

compList = c("CTRL-IDH")  
library(NetAct)
de.results = MicroDegs(eset = eset.brain_array) 

names(de.results)

# Save DE results
#----------------
save(de.results, file = paste('de.results.rda', sep = '', collapse = NULL))

# Load old gene set DB: hgs
#----------------
#fname.tfdb <- './hgs.rdata'  
#load(fname.tfdb)

#fname.tfdb <- './hgs.rcis.combined.filtered.rds' 
fname.tfdb <- './hgs.rds' 
hgs <- readRDS(file = fname.tfdb)
length(names(hgs))

#--------------------------
# Find TFs by GSEA analysis
#--------------------------
names(de.results)
tfs.CTRL.IDH = TF_GSEA(GSDB=hgs, #GSDB=hDB, 
              DErslt = de.results, #de.results$`CTRL-IDH`, #de.results$Overall, 
              minSize=8, 
              nperm = NO_PERMUTATIONS, 
              qval = T) 

names(tfs.CTRL.IDH)
class(tfs.CTRL.IDH)

# save TFs
fname.tfs <- './tfs.CTRL.IDH.csv'
write.csv(tfs.CTRL.IDH, file = fname.tfs, row.names = FALSE, quote = FALSE)
 
