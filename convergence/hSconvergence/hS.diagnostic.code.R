obtain.null.dist <- function(dataSimulation, dataReference, clusterCut, 
                             permutations, refClusterVar, corMethod, 
                             simulatedClusterVar){ 
  nClusters <- length(unique(clusterCut))
  nModelsKO <- dim(dataSimulation)[2] 
  n.gene <- dim(dataSimulation)[1]
  pValueMat <- matrix(0, nrow = nModelsKO, ncol = nClusters)  
  # create a matrix randomModels of size n.gene x permuations 
  # with the assigned values 1, 2, ..., n.gene repeated permutations 
  # times:  
  randomModels <-  matrix(rep(seq_len(n.gene), permutations), 
                          nrow = n.gene, ncol = permutations)  
  tmp <- randomModels
  # randomize the matrix randomModels:
  randomModels <- apply(randomModels, 2, sample) 

  # Calculate correlations between randomModels and dataReference: 
  permutedRefCor <- cor(randomModels, dataReference, method = corMethod)  
  
  # 
  for(j in seq_len(nClusters)){
    # Convert the correlations to distances: 
    # calculate the distances between each random model in randomModels and 
    # each of the samples (in dataRefence) belonging to the current cluster 
    # (square((1- cor(v1, v2))/2) is used as a distance 
    # between two vectors v1 and v2): 
    dist.mat <- ((1 - permutedRefCor[,which(clusterCut == j)])/2)^2 
    
    # (1) calculate the minimum distance between each sample and the random models
    # (2) sort these minimum distances 
    tempVector <- sort(apply(dist.mat, 1, .ClustFunction))
    for (i in seq_len(nModelsKO)){ 
      pValueMat[i,j] <- (which(abs(tempVector - simulatedClusterVar[i,j]) == min(abs(
        tempVector - simulatedClusterVar[i,j])))[1] - 1)/permutations 
      #[1] as sometimes which() might satisfy for multiple values 
    }
  }
  
  ret.obj <- list(permutedRefCor=permutedRefCor, pValueMat=pValueMat)
  
  return(ret.obj)
}


# v = obtain.null.dist(dataSimulation, dataReference, clusterCut, 
#                      permutations, refClusterVar, corMethod, 
#                      simulatedClusterVar)

# noRepeats <- 10
# permutations.v <- seq(1, noRepeats)*1000 

permutations.v <- c(1000, 5000, 10000, 15000, 20000, 40000, 60000, 80000, 100000)
noRepeats <- length(permutations.v)
retObj.list <- lapply(permutations.v, function(pv){
  return(obtain.null.dist(dataSimulation, dataReference, clusterCut, 
                          permutations=pv, refClusterVar, corMethod, 
                          simulatedClusterVar))
})



for(k in 1:noRepeats){
  #print(length(retObj.list[[k]]))
  print(dim(retObj.list[[k]]$permutedRefCor))
}

figdir <- paste('./figs-', circuit_idx, '/', sep = '')
dir.create(figdir) 


WIDTH <- 8  
HEIGHT <- 8 
figname <- paste(figdir, 'null.distributions-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
#par(mfrow = c(5, 2))   
par(mfrow = c(3, 3)) 
par(oma=c(3,3,3,3)) # b, l, t, r - all sides have 3 lines of space - outer margin
par(mar=c(1,1,4,1) + 0.1) # b, l, t, r - inner margin
xlimit <- c(-0.60, 0.6)
for(k in 1:noRepeats){
  #print(length(retObj.list[[k]])) 
  main.str <- as.character(permutations.v[k])
  hist(retObj.list[[k]]$permutedRefCor, breaks = 100, main = main.str, xlim = xlimit)
} 
dev.off()  



WIDTH <- 8  
HEIGHT <- 12
figname <- paste(figdir, 'pValues-CL1-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
#par(mfrow = c(5, 2))   
par(mfcol = c(9, 1)) 
par(oma=c(3,3,3,3)) # b, l, t, r - all sides have 3 lines of space - outer margin
par(mar=c(1,1,1,1) + 0.1) # b, l, t, r - inner margin

for(k in 1:noRepeats){
  #print(length(retObj.list[[k]])) 
  main.str <- as.character(permutations.v[k])
  plot(retObj.list[[k]]$pValueMat[,1], xlab='', xaxt='n', main=main.str)
} 
dev.off()  



WIDTH <- 12
HEIGHT <- 8
figname <- paste(figdir, 'pValues-CL1-1K-vs-100K-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
pdf(file = figname, width = WIDTH, height = HEIGHT, paper = 'special')
#par(mfrow = c(5, 2))   
par(mfcol = c(2, 1)) 
par(oma=c(3,3,3,3)) # b, l, t, r - all sides have 3 lines of space - outer margin
par(mar=c(1,1,1,1) + 0.1) # b, l, t, r - inner margin

main.str <- as.character(permutations.v[1])
plot(retObj.list[[1]]$pValueMat[,1], xlab='', xaxt='n', main=main.str)
main.str <- as.character(permutations.v[9])
plot(retObj.list[[9]]$pValueMat[,1], xlab='', xaxt='n', main=main.str)

dev.off()  

