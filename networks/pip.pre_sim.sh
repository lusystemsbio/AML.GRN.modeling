#!/usr/bin/env bash 
# Step 1. Create TF-target DBs - ../databases  

# Step 2. Generate core TFs by obtaining top TFs (5, 10, ..., 50) 
# from each of the three methods: (1) netact, (2) RI and (3) MARINa 
# and then combining them - ../tfSets/
# Output file: coreTFs.rds

date
# Step 3. Construct and simulate networks, and 
# create network metrics
R CMD BATCH inferCircuits.R
R CMD BATCH calCircuitSummary.R
R CMD BATCH selUniqCircuits.R

# Step 4. Calculate activities of TF in each circuit
R CMD BATCH calCircuit_TFactivities.R
date
