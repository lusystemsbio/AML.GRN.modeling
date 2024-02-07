

d= ((1 - permutedRefCor[,1])/2)^2 
length(d)


m=sort(apply(dist.mat, 1, min))
m
length(m)
head(m)
tail(m)

j=1 # untreated samples 
i=1 # first model
pval.tmp <- (which(abs(tempVector - simulatedClusterVar[i,j]) == min(abs(
  tempVector - simulatedClusterVar[i,j])))[1] - 1)/permutations

abs(tempVector - simulatedClusterVar[i,j])
min(abs(tempVector - simulatedClusterVar[i,j]))      

which(abs(tempVector - simulatedClusterVar[i,j]) == min(abs(tempVector - simulatedClusterVar[i,j])))[1]


simulated.cluster[which(simulated.cluster[,2] > pValue), 1] <- 0
