#!/usr/bin/env bash 

R CMD BATCH plotAccuracy.R
R CMD BATCH plotFlexibility.R
R CMD BATCH plot.rankings.flexibility-vs-accuracy.R
