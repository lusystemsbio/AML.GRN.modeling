#!/usr/bin/env bash 
date
R CMD BATCH xtract.sim.data.R
R CMD BATCH cal.mean_sd.R
R CMD BATCH unnormalize.data.R
R CMD BATCH build.model.R
date

# extra scripts
#R CMD BATCH plot.ACT_EXP.R
#R CMD BATCH plot.heatmap.ref_sim.R
#R CMD BATCH plotTopCircuitTPOs.R
#R CMD BATCH plot-sim.data-onPCs.R
