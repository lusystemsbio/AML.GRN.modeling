#!/usr/bin/env bash 

date 
# simulate circuits sequentially and save them as a list: 
#R CMD BATCH sim.circuits.R

# simulate circuits parallelly and save them as a list:
# (for 10K RACIPE models, performed on the server. 
# for 10K RACIPE models, the file size is too large to 
# handle and the subsequent script was written to save 
# the results separate file.)
#R CMD BATCH sim.circuits.para.R

# simulate circuits parallelly and save them as a separate file:
# (for 10K RACIPE models, performed on the server)
#R CMD BATCH sim.circuits.para.savedbyfile.R
date


