#!/usr/bin/env bash 

# Combine skd and dkd data sets:
# (1) Combine the proportion of models from skd and dkd 
#     across 100 repeats, sort them by the sorted mean of 
#     CL1 and save them
# (3) Plot bar plots for the proportions
# outputs: ./data/cluster_props.sorted.csv
R CMD BATCH comb.cluster_props.R

# For the combined data, perform Chi-squared test 
# inputs: ./data/cluster_props.sorted.csv
# outputs: ./data/stats.Chi-squared_test.csv
R CMD BATCH perform.chiSquared_test.R

# Plot the results from Chi-squared test
R CMD BATCH plot.stats.ChiSquared_test.R

# Calculate fold change for AML vs Untreated: ratio AML/Untreated
# inputs: ./data/cluster_props.sorted.csv
# outputs: ./data/foldChange.csv
R CMD BATCH cal.foldChange.R

# Plot fold change
# inputs: ./data/foldChange.csv
R CMD BATCH plot.foldChange.R
