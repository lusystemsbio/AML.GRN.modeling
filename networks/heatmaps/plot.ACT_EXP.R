#!/usr/bin/env Rscript

# Set working direcotry
#---------------------
rm(list=ls()) 
setwd(getwd())

# Set directories 
#---------------- 
figdir <- './figs/'
#dir.create(figdir)
 
# Load library functions
#----------------------- 
library(NetAct)
# library(Biobase) 
# library(dplyr)

# Load topology and obtain TFs
#----------------------------- 
fname_tpo  <- "./data/circuit-0.16-28-0.85.tpo"
tpo.df <- read.table(file = fname_tpo, header = T) 
tfs <- unique(union(as.character(tpo.df$Source), as.character(tpo.df$Target)))
length(tfs)

# Load TF activities
#-------------------
circuit_metrics <- read.csv(file = '../results/summary.circuits.sim.sorted.csv', row.names = 1)
circuit_idx <- rownames(circuit_metrics)[1]
circuit_idx # "0.16-28-0.85" - top circuit in all topologies

circuit_TFactivities.list <- readRDS(file = '../circuits/circuit_TFactivities.rds') 
names(circuit_TFactivities.list)
length(names(circuit_TFactivities.list))

circuit_TFactivities <- circuit_TFactivities.list[[circuit_idx]]
dim(circuit_TFactivities) # 73 x 20

# Retain only the ones that are in the final circuit
#---------------------------------------------------
circuit_TFactivities <- circuit_TFactivities[tfs, ]
dim(circuit_TFactivities) # 42 x 20

# Load expression file
#--------------------- 
fname.eset <- '../../data.tfs/eset.brain_array.rda'
load(fname.eset) 
class(eset.brain_array)

# Plot activities and expressions
#--------------------------------
WIDTH= 12 #10 #8 #12 #16 
HEIGHT=8 #10 #8 #12 #20 #12
fname_heatmap <- paste(figdir, 'activity_vs_expression.',  
                       toString(WIDTH), 'X', toString(HEIGHT), 
                       '.pdf', 
                       sep = '') 
pdf(file = fname_heatmap, width = WIDTH, height = HEIGHT, paper = 'special')
#Combine_heatmap(a$all_activities, eset=eset.brain_array) 
Combine_heatmap(circuit_TFactivities, eset=eset.brain_array) 
dev.off()

