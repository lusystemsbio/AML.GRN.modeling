
# sort by accuracy and flexibility
sortByTwoIndices.acc_and_flex <- function(circuit_metrics.sim){
  # circuit_metrics.sim <- circuit_metrics.sim.list[[sample_name]] 
  
  # sort by Accuracy:
  data.sortedByAccuracy <- circuit_metrics.sim[order(circuit_metrics.sim$Accuracy, decreasing = T), ]
  data.sortedByAccuracy$idxAccuracy <- 1:dim(data.sortedByAccuracy)[1]
  
  # # sort by AvgDist:
  # data.sortedByAvgDist <- data.sortedByAccuracy[order(data.sortedByAccuracy$AvgDist, decreasing = F), ]
  # data.sortedByAvgDist$idxAvgDist <- 1:dim(data.sortedByAvgDist)[1]
  
  # sort by flexibility:
  data.sortedByFlexibility <- data.sortedByAccuracy[order(data.sortedByAccuracy$flexibility, decreasing = T), ]
  data.sortedByFlexibility$idxFlexibility <- 1:dim(data.sortedByFlexibility)[1] # adding index
  
  # add both indices:
  data.sortedByFlexibility$idxBoth <- data.sortedByFlexibility$idxAccuracy+data.sortedByFlexibility$idxFlexibility
  
  # sort by conjugate index:
  data.sortedByBoth <- data.sortedByFlexibility[order(data.sortedByFlexibility$idxBoth, decreasing = F),]   
  return(data.sortedByBoth)
}










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

