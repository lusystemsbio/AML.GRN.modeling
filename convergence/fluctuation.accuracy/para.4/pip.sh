#!/usr/bin/env bash 

date
#R CMD BATCH cal.accuracy.para.R # run in the server
#R CMD BATCH draw.hist.accuracies.R
R CMD BATCH draw.hist.acc-3runs.R
date
