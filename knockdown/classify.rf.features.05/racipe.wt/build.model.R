#!/usr/bin/env Rscript

remove(list = ls()) 

# library(tidyverse)
#library(glmnet)
library(randomForest)

NO_CLUSTERS.REF <- 2 # 3 

outdir  <- './data/'
dir.create(outdir )

# Load training data 
#-------------------
NO_META_COLS = 1
fname_data <- './data/racipe.models.wt.labeled.csv'
mydata = read.csv(file=fname_data, header=TRUE)
dim(mydata)
mydata <- mydata[mydata$CLUSTER_NO!=3, ]
dim(mydata)

mydata$CLUSTER_NO <- as.numeric(mydata$CLUSTER_NO==1)
mydata$CLUSTER_NO <- factor(mydata$CLUSTER_NO)
unique(mydata$CLUSTER_NO)

# randomize data 
set.seed(1)
mydata <- mydata[sample(nrow(mydata)), ]

rf <- randomForest(CLUSTER_NO ~ MYC+TP53+RARA+RB1+ETS2,
                   data=mydata
                   )

fname.out <- paste(outdir, 'model.randomForest.rds', sep = '') 
saveRDS(rf, file = fname.out)



