#!/usr/bin/env bash 

R CMD BATCH remove.smallCircuits.R
R CMD BATCH replace.accuracies.R
R CMD BATCH rankSimCircuits.R
