#!/usr/bin/env bash 

date
# run on the server
#R CMD BATCH cal.metrics.sim_circuits.R
R CMD BATCH rankCircuits.R
R CMD BATCH cal.mean.sd.R
R CMD BATCH dotplot.top10circuits.R
date
