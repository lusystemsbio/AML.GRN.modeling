
normalize_by_wt_mean_and_sd <- function(geneExpression, means, sds){
   geneExpression <- sweep(geneExpression, 1, means, FUN = "-")
   geneExpression <- sweep(geneExpression, 1, sds, FUN = "/")  
   return(geneExpression)
}


obtain_clusterCut <- function(ndata, num_clusters){ 
   # ndata: rows are observations, columns are variables
   # corMethod <- 'pearson' # 'spearman' #
   # refCor <- cor((ndata), method = corMethod)
   # data.dist <- as.dist((1-refCor)/2)
   
   data.dist <- dist(t(ndata), method="euclidean")
   #data.dist <- dist(t(ndata), method="manhattan")
   #data.dist <- dist(t(ndata), method="minkowski")
   #data.dist <- dist(t(ndata), method="maximum")
   #data.dist <- dist(t(ndata), method="canberra")
   
   hc.rows <- hclust(data.dist, "ward.D2") #cluster along rows
   
   cutree_ranks = cutree(hc.rows, k=num_clusters) 
   return(cutree_ranks)
} 


obtain_data_by_cluster <- function(data.sim, clusterCut){
   data_by_cluster <- list() 
   for(cluster_no in 1:length(unique(clusterCut))){
      idx <- paste('CL', cluster_no, sep = '')
      data_by_cluster[[idx]] <- data.sim[, clusterCut==cluster_no]
   }   
   return(data_by_cluster)
}

