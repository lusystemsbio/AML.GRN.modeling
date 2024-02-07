#!/usr/bin/env Rscript

remove(list = ls()) 

# Common TFs among top 3 circuits (sorted by 4 way sorting)
#----------------------------------------------------------
tfs.ckt1 <- read.csv(file = './0.07-32-0.75/heatmaps/data/TFs-0.07-32-0.75.txt')
tfs.ckt2 <- read.csv(file = './0.07-32-0/heatmaps/data/TFs-0.07-32-0.txt')
tfs.ckt3 <- read.csv(file = './0.08-28-0.55/heatmaps/data/TFs-0.08-28-0.55.txt')

length(tfs.ckt1$x)
length(tfs.ckt2$x)
common.ckt1_ckt2 <- intersect(tfs.ckt1$x, tfs.ckt2$x)
length(common.ckt1_ckt2)

common.ckt1_ckt3 <- intersect(tfs.ckt1$x, tfs.ckt3$x)
length(common.ckt1_ckt3)

common.ckt2_ckt3 <- intersect(tfs.ckt2$x, tfs.ckt3$x)
length(common.ckt2_ckt3)

common.ckt1_ckt2_ckt3 <- intersect(common.ckt1_ckt2, tfs.ckt3$x)
length(common.ckt1_ckt2_ckt3)

# Common TFs among common TFs from above and TFs from flexible circuits
#----------------------------------------------------------------------
tfs.flex_ckt1 <- read.csv(file = './0.05-12-0/heatmaps/data/TFs-0.05-12-0.txt')
tfs.flex_ckt2 <- read.csv(file = './0.09-24-0.85/heatmaps/data/TFs-0.09-24-0.85.txt')

common.flex.ckt1_flex.ckt2 <- intersect(tfs.flex_ckt1$x, tfs.flex_ckt2$x)
length(common.flex.ckt1_flex.ckt2) # 14 

common.between.twoSets <-  intersect(common.ckt1_ckt2_ckt3, common.flex.ckt1_flex.ckt2)
length(common.between.twoSets) # 14

common.between.set1_and_flex.ctk1 <- intersect(common.ckt1_ckt2_ckt3, tfs.flex_ckt1$x) 
length(common.between.set1_and_flex.ctk1) # 20 

common.between.set1_and_flex.ctk2 <- intersect(common.ckt1_ckt2_ckt3, tfs.flex_ckt2$x) 
length(common.between.set1_and_flex.ctk2) # 22 


# Common TFs between top circuit in phase.46 and phase.37
#---------------------------------------------------------
tf.df <- read.csv('/Users/kateba/research/aml.idh/phase.37/networks/pathway.anno/tfs.annotated/TFsWith_singleAnnotedPathway.csv')

tfs.phase.37 <- tf.df$tf
class(tfs.phase.37) 
length(tfs.phase.37) # 42
length(tfs.ckt1$x) # 52
common.between.phase.37_and_46 <- intersect(tfs.phase.37, tfs.ckt1$x)
length(common.between.phase.37_and_46) # 40 


common.between.phase.37_and_46
setdiff(tfs.phase.37, tfs.ckt1$x)
sort(setdiff(tfs.ckt1$x, tfs.phase.37))


# Common TFs among top circuits (Acc, flex) vs (Acc, avg dist, flex)
#===================================================================================
tfs.ckt2.Acc_flex <- read.csv(file = './0.09-32-0.85/heatmaps/data/TFs-0.09-32-0.85.txt') 
length(tfs.ckt2.Acc_flex$x) # 29

tfs.ckt3.Acc_flex <- read.csv(file = './0.07-36-0.85/heatmaps/data/TFs-0.07-36-0.85.txt') 
length(tfs.ckt3.Acc_flex$x) # 32

# 2nd and 3rd top circuits (Acc, flex)
common.tfs <- intersect(tfs.ckt2.Acc_flex$x, tfs.ckt3.Acc_flex$x)
length(common.tfs)
setdiff(tfs.ckt2.Acc_flex$x, tfs.ckt3.Acc_flex$x)
setdiff(tfs.ckt3.Acc_flex$x, tfs.ckt2.Acc_flex$x)


# 2nd top circuit (Acc, flex) and 1st top circuit (Acc, avg dist, flex)
common.tfs <- intersect(tfs.ckt2.Acc_flex$x, tfs.ckt1$x)
length(common.tfs) # 28 
setdiff(tfs.ckt2.Acc_flex$x, tfs.ckt1$x)
length(setdiff(tfs.ckt1$x, tfs.ckt2.Acc_flex$x))


# 2nd top ckt (Acc, flex) and top ckt (phase.37)
common.tfs <- intersect(tfs.phase.37, tfs.ckt2.Acc_flex$x)
length(common.tfs)
length(setdiff(tfs.phase.37, tfs.ckt2.Acc_flex$x))
length(setdiff(tfs.ckt2.Acc_flex$x, tfs.phase.37))

