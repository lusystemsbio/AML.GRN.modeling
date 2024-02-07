#!/usr/bin/env bash 

date
R CMD BATCH comb.cluster_props.R
R CMD BATCH plot.cluster.props.R
date
