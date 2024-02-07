#!/usr/bin/env bash 

date
#R CMD BATCH cal.accuracy.para.R # run in the server
R CMD BATCH draw.hist.accuracies.R
date
