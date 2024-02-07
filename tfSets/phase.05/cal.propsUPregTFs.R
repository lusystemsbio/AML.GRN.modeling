# Input files: 
# (1) target DB - used for constructing initial network 
# (2) core TFs - used for extracting a possible topology from the target DB 
# (3) Expression data - raw_brainarray.sele.txt - used for 
#                       restricting the possible topology specific to data
 
remove(list = ls()) 

# Constants
#-------------
COLNAME.PREFIX.CONTROL <- 'CTRL' 
COLNAME.PREFIX.TREATMENT <- 'IDH'

library(Biobase) 
library(dplyr)
#library(NetAct)
source('./functions.tfBalance.R')

# Output directories
#-------------------
outdir <- './data/'
dir.create(outdir)

# Load brain array expression data as an ExpressionSet object
#------------------------------------------------------------
fname.eset.brain_array <- '../../data.tfs/eset.brain_array.rda'
load(fname.eset.brain_array) # loads object eset.brain_array 
class(eset.brain_array) 
data <- exprs(eset.brain_array) 
dim(data)

# Load DE results 
#---------------- 
fname.de.results <- '../../data.tfs/de.results.rda'
load(fname.de.results) # loads object de.results 


# select TOP n TFs from each of the THREE methods: NetAct, MARINa, and RI
#------------------------------------------------------------------------
coreTFs.list <- readRDS('./data/coreTFs.rds')
names(coreTFs.list)

# Construct possible circuit set for the TFs based on each TF-target DB in the TF-target DB set
#---------------------------------------------------------------------------------------------- 
targetDB.list <- readRDS('../../databases/targetDB.list.rds') 


# Calculate the number of UP TFs in CTRL and TRTMT
#---------------------------------------------------------------------------------
cNames <- c('feature_ratio', 'no_tfs', 'no_tfs_actual' , 'method', 'CTRL', 'TRTMT') 
# stats.upTF.df <- as.data.frame(matrix(nrow = 0, ncol = length(cNames)), stringsAsFactors = FALSE)

stats.upTF.df <- as.data.frame(matrix(nrow = (length(coreTFs.list)*length(coreTFs.list$`0.05`)*length(coreTFs.list$`0.05`$`5`)),
                                      ncol = length(cNames)),
                               stringsAsFactors = FALSE) 
colnames(stats.upTF.df) <- cNames

count <- 1 
for(fr in names(targetDB.list)){  
   cat('feature ratio cut off:')
   print(fr)
   targetDB.cur <- targetDB.list[[fr]] 
   coreTFs.byTOPtfCount <- coreTFs.list[[fr]]
   for(top.TFs.count in names(coreTFs.byTOPtfCount)){  
      coreTFs.byMethod <- coreTFs.byTOPtfCount[[top.TFs.count]] 
      
      for(method.name in names(coreTFs.byMethod)){ 
         tfs.sele <- coreTFs.byMethod[[method.name]] 
         upGeneCount <- cal.UPgeneCount(tfs = tfs.sele, 
                                        GSDB = targetDB.cur, #targetDB.list[[fr]], #targetDB, #hgs,
                                        eset = data,
                                        de.results = de.results,
                                        COLNAME.PREFIX.CONTROL = COLNAME.PREFIX.CONTROL,
                                        COLNAME.PREFIX.TREATMENT = COLNAME.PREFIX.TREATMENT
         )  
         stats.upTF.df[count, ] <- c(fr, top.TFs.count, length(tfs.sele), method.name, 
                                     upGeneCount[["CTRL"]], upGeneCount[["TRTMT"]])
         count <- count+1
         #break()
      }
      #break()
   }
   #break()
}

# Calculate the ratio of UP TFs between CTRL and TRTMT
ratio_upTFs <- as.numeric(stats.upTF.df$TRTMT)/as.numeric(stats.upTF.df$CTRL)
min(ratio_upTFs)
max(ratio_upTFs)
tmp.df <- cbind(stats.upTF.df, ratio_upTFs)
colnames(tmp.df) <- c(colnames(stats.upTF.df), 'TRTMT_CTRL')
stats.upTF.df <- tmp.df


fname.out <- paste(outdir, 'stats.upTF.csv', sep = '')
write.csv(stats.upTF.df, file = fname.out, quote = F, row.names = F)

method <- "Netact"
stats.df <- stats.upTF.df[stats.upTF.df$method==method,] 
fname.out <- paste(outdir, 'stats.', method ,'.csv', sep = '')
write.csv(stats.df, file = fname.out, quote = F, row.names = F) 

method <- "MARINa"
stats.df <- stats.upTF.df[stats.upTF.df$method==method,] 
fname.out <- paste(outdir, 'stats.', method ,'.csv', sep = '')
write.csv(stats.df, file = fname.out, quote = F, row.names = F) 


method <- "RI"
stats.df <- stats.upTF.df[stats.upTF.df$method==method,] 
fname.out <- paste(outdir, 'stats.', method ,'.csv', sep = '')
write.csv(stats.df, file = fname.out, quote = F, row.names = F) 

method <- "COMB"
stats.df <- stats.upTF.df[stats.upTF.df$method==method,] 
fname.out <- paste(outdir, 'stats.', method ,'.csv', sep = '')
write.csv(stats.df, file = fname.out, quote = F, row.names = F) 

