#!/usr/bin/env bash 
date
#R CMD BATCH gen.hS.object.R
R CMD BATCH cal.mean_sd.R
R CMD BATCH xtract.sim.data.R

#R CMD BATCH plot.ACT_EXP.R
#R CMD BATCH plot.heatmap.ref_sim.R
#R CMD BATCH plotTopCircuitTPOs.R
#R CMD BATCH plot-sim.data-onPCs.R
date
