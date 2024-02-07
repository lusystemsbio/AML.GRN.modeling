#!/usr/bin/env bash 
date 

# calculate avg distance for each circuit
# (using heatmap similarity)
# need to check the use of this function call???
# output: ./results/summary.circuits.sim.csv
R CMD BATCH cal.metrics.sim_circuits.R

# calculate mean and sd of avg distance for each circuit
# (calculation performed in the server)
R CMD BATCH cal.metrics.mean.sd.R

# calculate circuit flexibility: 
# output: circuit.flexibility
R CMD BATCH cal.circuitFlexibility.R

# attach calculated circuit flexibility at the end of summary.circuits.sim.csv:
# summary.circuits.sim.flex.csv
R CMD BATCH attach.circuitFlexibility.R

# attach significance value to each circuit:
# input:./results/summary.circuits.sim.flex.csv
# output: ./summary.circuits.sim.sigVal.csv
R CMD BATCH cal.maxSigValByMethod.R

# rank circuits based on 3way sorting: 
# input: ./summary.circuits.sim.sigVal.csv
# output: ./results/summary.circuits.sim.sorted.csv
#R CMD BATCH rankSimCircuits.R

# find how many TFs from each method are included from each method
# input: ./results/summary.circuits.sim.sorted.csv
# output: ./results/tfCount.byMetod.csv
#R CMD BATCH tfCountByMethod.R

R CMD BATCH rankSimCircuits.R

# remove small circuits (fewer than 15 nodes) from file sorted by ACC and FLEX
R CMD BATCH remove.small.circuits.R
date
