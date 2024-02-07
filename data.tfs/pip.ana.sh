#!/usr/bin/env bash

date 
R CMD BATCH evalTOP.TFs.R
R CMD BATCH cal.ACT_EXP.R
date 

