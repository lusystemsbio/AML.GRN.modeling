
# Given a set of genes, what are the TFs they are targetted by
#------------------------------------------------------------
mapGenes_to_TFs <- function(targetDB, sigGenes){
   TFs.mapped <- c()
   for(tf in names(targetDB)){
      targets <- targetDB[[tf]] 
      commonGenes <- intersect(sigGenes, targets) 
      if(length(commonGenes) != 0){
         TFs.mapped <- c(TFs.mapped, tf)
      }
   }  
   return(TFs.mapped)
}

# Fisher Exact Test
#-------------------------
FET_1not2 = function(glist1, glist2, ntot){
   a = sum(sign(match(glist1, glist2, nomatch = 0)))
   b = length(glist2) - a
   c = length(glist1) - a   #actual counts
   d = ntot - a - b - c
   tmp_mat = matrix(c(c, d, a, b), ncol = 2)
   tmp_fisher = fisher.test(tmp_mat)
   tmp_pval = tmp_fisher$p.value
   exp = (a+c)*(c+d)/ntot
   #  cat(exp,c, "\n")
   return(list(pval = tmp_pval, over_expect = ifelse(c > exp, 1, 0), act = c, exp = exp))
}


# Calculate average for sub samples of simulated data
cal.avgDist.bySample <- function(data.REF, data.sim, clusterCut.REF, SAMPLE_SIZE){
   NO_SAMPLES <- dim(data.sim)[2]/SAMPLE_SIZE
   
   AvgDist <- vector(mode = 'numeric', length = NO_SAMPLES)
   for(SAMPLE_NO in 1:NO_SAMPLES){
      start.idx <- SAMPLE_SIZE * (SAMPLE_NO-1) + 1
      end.idx <- SAMPLE_SIZE * SAMPLE_NO
      # print(start.idx)
      # print(end.idx)
      data.sim.sub <- data.sim[, start.idx:end.idx] 
      # print(dim(data.sim.sub)) 
      hS <- sracipeHeatmapSimilarity(dataReference = data.REF,
                                     #dataSimulation = data.sim, 
                                     dataSimulation = data.sim.sub,  
                                     returnData = F, 
                                     #nClusters = NUM_CLUSTERS, 
                                     clusterCut = clusterCut.REF)
      AvgDist[SAMPLE_NO] <- hS$AvgDist
   } 
   return(AvgDist)
}

