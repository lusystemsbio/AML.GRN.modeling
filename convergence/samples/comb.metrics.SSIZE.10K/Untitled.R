perc.rank <- function(x) trunc(rank(x))/length(x)


my.df <- data.frame(x=rnorm(200))
my.df <- within(my.df, xr <- perc.rank(x))
