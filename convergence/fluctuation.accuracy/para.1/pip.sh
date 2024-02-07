#!/usr/bin/env bash 

date
#R CMD BATCH cal.accuracy.seq.R
#R CMD BATCH cal.accuracy.para.R
R CMD BATCH draw.hist.accuracies.R
date
