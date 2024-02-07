
# Set working direcotry
#---------------------
rm(list=ls()) 
setwd(getwd())

NUM_OF_TFS <- 50

# Set directories 
#---------------- 
libdir <- '../lib/'

figdir <- './figs/'
dir.create(figdir)
 
# Load library functions
#----------------------- 
library(NetAct) 
library(Biobase) 
library(dplyr)

# Load TFs file
#-------------- 
fname_tf  <- "./tfs.Netact.txt"
tfs.all.df <- read.table(file = fname_tf, header = T)   

tfs.df <- data.frame(tfs.all.df$tf[1:NUM_OF_TFS]) 
colnames(tfs.df) <- 'TF'

# Load expression file
#--------------------- 
fname.eset <- 'eset.brain_array.rda'
load(fname.eset) 
class(eset.brain_array) 

# Load DE results
#--------------------- 
fname.de.results <- 'de.results.rda'
load(fname.de.results) 

# Obtain expressions for CONTROL and TREATMENT 
#--------------------------------------------
# Get data for the current comparison 
data <- exprs(eset.brain_array) 
dim(data) 

# Load old gene set DB: hgs
#----------------
# fname.tfdb <- './hgs.rdata'
# load(fname.tfdb)
# length(names(hgs))

#fname.tfdb <- './hgs.RI.combined.filtered.rds' 
#fname.tfdb <- './hgs.rcis.combined.filtered.rds' 
fname.tfdb <- './hgs.rds' 
hgs <- readRDS(fname.tfdb)
length(names(hgs))

names(de.results)

# Retain the TFs that are in the TF-target DB - used when 
#-------------------------------------------
tfs.valid <- intersect(names(hgs), tfs.df$TF)
length(tfs.valid)  

# Calculate activities 
#---------------------- 
# a = TF_Activity(tfs = as.character(tfs.df$TF),
#                 GSDB=hgs, #GSDB=hDB,
#                 eset=data,
#                 DErslt=de.results$Overall #de.results  #de.results$Overall
# )

a = TF_Activity(tfs = as.character(tfs.df$TF), #tfs.valid 
                GSDB=hgs, #GSDB=hDB,
                eset=data,
                DErslt=de.results  #de.results$Overall
)

class(a)
dim(a$all_activities)
act.tmp <- a$all_activities
 
# Plot activities and expressions
#--------------------------------
WIDTH= 12 #16 
HEIGHT=12 #20 #12
fname_heatmap <- paste(figdir, 'activity_vs_expression.',  
                       toString(WIDTH), 'X', toString(HEIGHT), 
                       '.pdf', 
                       sep = '') 
pdf(file = fname_heatmap, width = WIDTH, height = HEIGHT, paper = 'special')
Combine_heatmap(a$all_activities, eset=eset.brain_array) 
dev.off()

