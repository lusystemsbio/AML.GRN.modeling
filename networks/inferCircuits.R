# Input files: 
# (1) target DB - used for constructing initial network 
# (2) core TFs - used for extracting a possible topology from the target DB 
# (3) Expression data - raw_brainarray.sele.txt - used for 
#                       restricting the possible topology specific to data
 
# Output file(s): 
# Constructed network topologies in two formats  
#  (1) circuits: 
#      Each entry is a circuit   
#      Index of the list:  
#        a conjugate key in the format: TopTFs-FeatureRatio-AbsCor 
#      There is one circuit for each conjugate key


remove(list = ls()) 

# List of interaction strengths
#INTERACTION.STRENGTHS  <- seq(0.05, 0.95, 0.05)
INTERACTION.STRENGTHS  <- seq(0.00, 0.95, 0.05)

# Internal parameters
SUBNETWORK.SIZE.TSH <- 0.80 


library(igraph)
source('./functions.R')
library(sRACIPE)
library(NetAct)


# Output directories
#-------------------
outdir <- './circuits/'
dir.create(outdir)


# Load brain array expression data as an ExpressionSet object
#------------------------------------------------------------
fname.eset.brain_array <- '../data.tfs/eset.brain_array.rda'
load(fname.eset.brain_array) # loads object eset.brain_array 


# Load DE results 
#---------------- 
fname.de.results <- '../data.tfs/de.results.rda'
load(fname.de.results) # loads object de.results 


# select TOP n TFs from each of the THREE methods: NetAct, MARINa, and RI
#------------------------------------------------------------------------
coreTFs.list <- readRDS('../tfSets/data/coreTFs.rds')
names(coreTFs.list)

# Load TF-target DBs
#-------------------  
targetDB.list <- readRDS('../databases/targetDB.list.rds')


# Infer Circuits
#----------------
circuits <- list()   
for(fr in names(coreTFs.list)){ 
   print(fr)
   coreTFs.byTOPtfCount <- coreTFs.list[[fr]] 
   for(top.TFs.count in names(coreTFs.byTOPtfCount)){
      print(top.TFs.count) 
      coreTFs.byMethod <- coreTFs.byTOPtfCount[[top.TFs.count]] 
      coreTFs <- coreTFs.byMethod[["COMB"]] # extract combined TFs 
      
      circuits.tmp <- infer.circuits.by.targetDB(coreTFs = coreTFs,
                                                 targetDB = targetDB.list[[fr]], #targetDB,
                                                 eset.brain_array = eset.brain_array,
                                                 de.results = de.results, 
                                                 int.strengths <- INTERACTION.STRENGTHS, 
                                                 SUBNETWORK.SIZE.TSH=SUBNETWORK.SIZE.TSH)  
      # extract and save individual circuits: 
      for(mi.tsh in names(circuits.tmp)){
         idx_name <- paste(fr, top.TFs.count, mi.tsh, sep = '-')
         print(idx_name)
         circuits[[idx_name]] <- circuits.tmp[[mi.tsh]]
      }
      #break()
   } 
   #break()
}

names(circuits)
length(names(circuits))
head(names(circuits))

fname.circuits <- paste0(outdir, 'circuits.rds')
saveRDS(circuits, file = fname.circuits) 