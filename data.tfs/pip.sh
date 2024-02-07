#!/usr/bin/env bash
date 

# obtain TFs using Netact method: 
R CMD BATCH find_TFs.R
R CMD BATCH order.TFs.R


# Combine TOP TFs from three method-Netact, MARINa, and RI: 
R CMD BATCH comb.TFs.from_three_methods.R
date 