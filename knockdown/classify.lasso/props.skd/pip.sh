#!/usr/bin/env bash 

R CMD BATCH cal.cluster.props.R
R CMD BATCH plot.cluster.props.R
R CMD BATCH plot.kd_clusters_on_PCs.R