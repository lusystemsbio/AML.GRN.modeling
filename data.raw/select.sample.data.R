"
 (1) Load raw brain array data 
 (2) Load and combine ids for control and two treatments  
 (3) Select raw brain array data for the combined ids 
 (4) Change columans of the selected data with treatment tag
 (4) Save selected raw brain array data
"

rm(list=ls())
setwd(getwd())

#------------------------------------------------------#
# CONSTANTS
#------------------------------------------------------#
# TREATMENT.1.NAME <- 'IDH1'
# TREATMENT.2.NAME <- 'IDH2'

TREATMENT.1.NAME <- 'IDH'
TREATMENT.2.NAME <- 'IDH'

TREATMENT.NAME <- 'IDH' # includes IDH1 and IDH2
#------------------------------------------------------#
# Load libray functions
#------------------------------------------------------#
# libdir <- '../lib/'
source('./lib.data.raw.R')

#------------------------------------------------------#
# Load and map the column names of raw AML data
#------------------------------------------------------#
# INPUT FILE NAMES: 
# raw AML expression file:
fname_data <- "raw_brainarray.txt" 
# sample types: control vs treatment
fname_mapping ="erasmus_exp.txt"
rdata <- map.colnames_of_raw_data(fname_data = fname_data, 
                                  fname_mapping = fname_mapping)

#------------------------------------------------------#
# Load and combine CONTROL and TREATMENT groups
#------------------------------------------------------#

# File names for CONTROL AND TREATMENT groups
#-------------------------------------------#
fname_samples.CTRL <- 'samples.CONTROL.csv'  
fname_samples.TREATMENT.1 <- 'IDH1.txt'  
fname_samples.TREATMENT.2 <- 'IDH2.txt' 


# Load control ids:
samples.CTRL <- read.csv(fname_samples.CTRL, row.names = 1)

# Load ids treatment 1: 
samples.TREATMENT.1 <- read.table(file = fname_samples.TREATMENT.1, 
                                  col.names = c("PATIENT.ID")) 

# Load ids treatment 2: 
samples.TREATMENT.2 <- read.table(file = fname_samples.TREATMENT.2, 
                                  col.names = c("PATIENT.ID")) 

#samples.TREATMENT.1 <- as.character(samples.TREATMENT.1$PATIENT.ID)


sample_ids.TREATMENT.1 <- samples.TREATMENT.1$PATIENT.ID  
sample_ids.TREATMENT.2 <- samples.TREATMENT.2$PATIENT.ID 

# Load treatment IDs:
#--------------------
fname_samples.TREATMENT <- './IDH.txt'  
samples.TREATMENT <- read.table(file = fname_samples.TREATMENT, 
                                col.names = c("PATIENT.ID"))
sample_ids.TREATMENT <- samples.TREATMENT$PATIENT.ID  

# create column names and save them
#---------------------------------
colnames.new <- create.colnames(samples.CTRL, 
                                sample_ids.TREATMENT,
                                TREATMENT.NAME) 
fname_samples <- 'samples.selected.txt'
write.table(colnames.new,file = fname_samples, quote = FALSE, 
            row.names = FALSE, col.names = c("SAMPLE.ID")) 

# Combine sample names:
#-----------------------#
sample_ids.comb <- bind.sample_names(samples.CTRL, 
                                     sample_ids.TREATMENT)

length(sample_ids.comb)
sample_ids.comb


# Extract data for the selected samples
#----------------------------------------#
rdata.sele <- rdata[ ,sample_ids.comb]
dim(rdata.sele)

# Change column names of the selected data with column tag 
# of TREATMENT TYPE/PATIENT TYPE
#-----------------------------------------------------------#
names(rdata.sele) <- colnames.new


# Save the selected data
#-------------------------#
cat('Save selected brain array data')
fname_data.sele <- 'raw_brainarray.sele.txt'
write.table(rdata.sele, file = fname_data.sele, sep='\t',
            quote = FALSE, row.names = TRUE)
