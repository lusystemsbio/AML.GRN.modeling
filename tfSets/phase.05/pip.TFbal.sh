#!/usr/bin/env bash 
date
R CMD BATCH selTopTFs.byMethod.R
R CMD BATCH cal.propsUPregTFs.R
R CMD BATCH plotBalTFs.R
R CMD BATCH dotplot.RatioUPTFs.R
date