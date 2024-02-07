
figdir <- './figs.ri/'
dir.create(figdir)
 
# Load library functions
#----------------------- 
library(NetAct) 
library(Biobase) 
library(dplyr)

# Load TFs file
#-------------- 
# obtain  TFs.xwen from xplore.TFs.RI.R 


# Load expression file
#--------------------- 
fname.eset <- '../../data.tfs/eset.brain_array.rda'
load(fname.eset) 
class(eset.brain_array) 

# Load DE results
#---------------- 
fname.de.results <- '../../data.tfs/de.results.rda'
load(fname.de.results) 

# Obtain expressions for CONTROL and TREATMENT 
#--------------------------------------------
# Get data for the current comparison 
data <- exprs(eset.brain_array) 
dim(data) 

# Load old gene set DB: hgs
#--------------------------
targetDB.list <- readRDS('../targetDB.list.rds')
hgs <- targetDB.list$`0.05`
names(de.results)

# Retain the TFs that are in the TF-target DB - used when 
#--------------------------------------------------------
tfs.valid <- intersect(names(hgs), TFs.xwen)
length(tfs.valid)  

# Calculate activities 
#----------------------
a = TF_Activity(tfs = tfs.valid,
                GSDB=hgs, #GSDB=hDB,
                eset=data,
                DErslt=de.results  #de.results$Overall
)

class(a)
dim(a$all_activities)
act.tmp <- a$all_activities

# Plot activities and expressions
#--------------------------------
WIDTH=12 #16 
HEIGHT=6 #12 #20 #12
fname_heatmap <- paste(figdir, 'activity_vs_expression.',  
                       toString(WIDTH), 'X', toString(HEIGHT), 
                       '.pdf', 
                       sep = '') 
pdf(file = fname_heatmap, width = WIDTH, height = HEIGHT, paper = 'special')
Combine_heatmap(a$all_activities, eset=eset.brain_array) 
dev.off()

