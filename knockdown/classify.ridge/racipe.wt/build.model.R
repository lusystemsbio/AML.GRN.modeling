#!/usr/bin/env Rscript

remove(list = ls()) 

# library(tidyverse)
library(glmnet)

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

x <- model.matrix(CLUSTER_NO ~., mydata)[,-1] 
y <- ifelse(mydata$CLUSTER_NO == 1, 1, 0)
head(y)


ALPHA <- 0 # Ridge regression

# perform 10 fold cross validation to obtain best lamba (penalty) 
cv.glmnet.obj <- cv.glmnet(x, y, alpha = ALPHA, family = "binomial") 

# create glmnet model using the optimum lambda
model.glmnet <- glmnet(x, y, alpha = ALPHA, family = "binomial",
                lambda = cv.glmnet.obj$lambda.1se #cv.glmnet.obj$lambda.min
)

fname.out <- paste(outdir, 'cv.glmnet.obj.rds', sep = '') 
saveRDS(cv.glmnet.obj, file = fname.out)

fname.out <- paste(outdir, 'model.glmnet.rds', sep = '') 
saveRDS(model.glmnet, file = fname.out)

