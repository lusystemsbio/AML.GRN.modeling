#!/usr/bin/env bash 

date
R CMD BATCH remove.smallCircuits.R
R CMD BATCH rankSimCircuits.R
date
