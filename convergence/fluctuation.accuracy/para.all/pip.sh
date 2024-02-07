#!/usr/bin/env bash 

date
#R CMD BATCH cal.accuracy.para.R # run in the server
R CMD BATCH draw.hist.accuracies.R
R CMD BATCH draw.hist.acc-4runs.R
R CMD BATCH cal.mean.sd.R
date
