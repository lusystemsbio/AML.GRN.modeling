#/usr/bin/env bash 

date 
#run on server: 
# calculate metrics: Accuracy and AvgDist:
R CMD BATCH cal.metrics.sim_circuits.R

# Rank the circuits: 3-way sorting
R CMD BATCH rankCircuits.R

# Calculate mean and sd of combined scores (ranks) across samples 
# for each circuit: 
R CMD BATCH cal.mean.sd.R
R CMD BATCH dotplot.sSize-10K.top10circuits.R
date 
