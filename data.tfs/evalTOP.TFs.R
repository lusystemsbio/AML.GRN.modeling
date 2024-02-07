#!/usr/bin/env Rscript
rm(list=ls()) 
setwd(getwd())

NUM_OF_TFS <- 50

# Load library functions
#----------------------- 
library(NetAct) 
library(Biobase) 
library(dplyr)

# Load TFs file
#-------------- 
fname_tf  <- "./tfs.Netact.txt"
tfs.all.df <- read.table(file = fname_tf, header = T) 

tfs.netact <- data.frame(tfs.all.df[1:NUM_OF_TFS,] ) 
stats.tf <- data.frame(as.character(tfs.netact$tf), tfs.netact$qvals, -log10(tfs.netact$qvals)) 
colnames(stats.tf) <- c('tf', 'qval', 'nlogqval')
class(stats.tf$qval)
class(stats.tf$tf)

library(ggplot2)

class(stats.tf$tf)
stats.tf$tf <- factor(stats.tf$tf, levels =  stats.tf$tf) 


# Calculate TF activities 
#------------------------
fname.tfdb <- './hgs.rds' 
hgs <- readRDS(fname.tfdb)
length(names(hgs))

# Load DE results
#--------------------- 
fname.de.results <- 'de.results.rda'
load(fname.de.results)

# Load expression file
#--------------------- 
fname.eset <- 'eset.brain_array.rda'
load(fname.eset) 
class(eset.brain_array) 

# Obtain expressions for CONTROL and TREATMENT 
#--------------------------------------------
# Get data for the current comparison 
data <- exprs(eset.brain_array) 
dim(data) 


# Retain the TFs that are in the TF-target DB
#-------------------------------------------
tfs.valid <- intersect(names(hgs), as.character(stats.tf$tf))
length(tfs.valid)

# Calculate activities 
#---------------------- 
a = TF_Activity(tfs = tfs.valid, # as.character(tfs.df$TF),
                GSDB=hgs, #GSDB=hDB,
                eset=data,
                DErslt=de.results  #de.results$Overall
)

tf.act <- a$all_activities

# check and remove lost TFs from stats.df 
tfs_to_remove <- setdiff(as.character(stats.tf$tf), rownames(tf.act))
tmp.df <- stats.tf[which(!(as.character(stats.tf$tf) %in% tfs_to_remove)), ]
dim(tmp.df)
stats.tf <- tmp.df

dim(tf.act)
tf.act.ordered <- tf.act[as.character(stats.tf$tf), ] 


tf.act_avg <- data.frame(rowMeans(tf.act.ordered[,1:11]),  rowMeans(tf.act.ordered[,12:20])) 
colnames(tf.act_avg) <- c('ctrl', 'aml')

tmp.df <- stats.tf
#stats.tf <- data.frame(tmp.df$tf, tmp.df$qval, tmp.df$nlogqval, tf.act_avg$ctrl, tf.act_avg$aml)
stats.tf <- data.frame(tmp.df, tf.act_avg$ctrl, tf.act_avg$aml)
colnames(stats.tf) <- c(colnames(tmp.df), 'act.ctrl', 'act.aml')


# Plot for q values
#-----------------
p1 <- ggplot(stats.tf) + 
   geom_bar(aes(x=tf, y=nlogqval),
            stat="identity", fill="forestgreen", alpha=0.5)
p1
p2 <- p1 + theme(axis.text.x = element_blank(),  #element_text(angle = 90, hjust = 0.7, vjust = 0.5),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank())
p2
plot.qval <- p2

# Plot for activity (CTRL)
#-------------------------
p1 <- ggplot(stats.tf) + 
   geom_bar(aes(x=tf, y=act.ctrl),
            stat="identity", fill="blue", alpha=0.5) + ylim(c(-0.90, 0.90))
p1
p2 <- p1 + theme(axis.text.x = element_blank(), #element_text(angle = 90, hjust = 0.7, vjust = 0.5), 
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank()) 

p2

plot.act.ctrl <- p2


# Plot for activity (AML)
#-------------------------
p1 <- ggplot(stats.tf) + 
   geom_bar(aes(x=tf, y=act.aml),
            stat="identity", fill="red", alpha=0.5) + ylim(c(-0.9, 0.9))
p1
p2 <- p1 + theme(axis.text.x = element_text(angle = 90, hjust = 0.7, vjust = 0.5), 
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank()) 
p2

plot.act.aml <- p2

plot.list <- list(plot.qval, plot.act.ctrl, plot.act.aml)
library(gridExtra)
plot_on_grid <- grid.arrange(grobs=plot.list, nrow=3)



WIDTH <- 12
HEIGHT <- 6
outdir <- './figs/'
dir.create(outdir)
fname.out <- paste(outdir, 'qvals-vs-activities-', WIDTH, '-', HEIGHT, '.pdf', sep = '')
ggsave(filename = fname.out, plot_on_grid, width = WIDTH, height = HEIGHT) 
