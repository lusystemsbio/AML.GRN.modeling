#!/usr/bin/env Rscript
remove(list = ls())


perc.rank <- function(x, xo)  length(x[x <= xo])/length(x)*100

perc.rank <- function(x) trunc(rank(x))/length(x)

my.df <- data.frame(x=rnorm(200))
my.df <- within(my.df, xr <- perc.rank(x))
my.df


