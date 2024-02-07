
# Calculate euclidan distance between each distribution of KD simulations and 
# WT distribution and return their average
calDistance <- function(racipe.kd){
  prop.wt <- racipe.kd$WT
  dist.df <- as.data.frame(matrix(nrow = (length(racipe.kd)-1), ncol = 2)) 
  colnames(dist.df) <- c('tf', 'dist')
  rownames(dist.df) <- names(racipe.kd)[2:length(racipe.kd)]
  for(tf in names(racipe.kd)[2:length(racipe.kd)]){
    #print(tf)
    prop.kd <- racipe.kd[[tf]] 
    d <- sqrt((as.numeric(prop.wt[1])-as.numeric(prop.kd[1]))^2 + 
                (as.numeric(prop.wt[2])-as.numeric(prop.kd[2]))^2)
    dist.df[tf,] <- c(tf, d)
    #break()
  }
  avg.dist <- sum(as.numeric((dist.df$dist)))/length(dist.df$dist)   
  return(avg.dist)
}
