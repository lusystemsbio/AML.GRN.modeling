normalize_by_wt_mean_and_sd <- function(geneExpression, means, sds){
   geneExpression <- sweep(geneExpression, 1, means, FUN = "-")
   geneExpression <- sweep(geneExpression, 1, sds, FUN = "/")  
   return(geneExpression)
}

cal.model_count_by_cluster <- function(predicted_cluster_no, NO_CLUSTERS.REF){
   model_count_by.cluster <- vector(mode = 'numeric', length = (length(NO_CLUSTERS.REF)+1))
   for(cl_no in 0:NO_CLUSTERS.REF){
      #print(cl_no) 
      model_count_by.cluster[cl_no+1] <- sum(predicted_cluster_no==cl_no)
   }
   return(model_count_by.cluster)
}

predict_cluster_no <- function(mymodel, ndata.kd, PROB_DIFF_THRESHOLD){
   predicted_ranks_l <- predict(mymodel, ndata.kd) #predicted_ranks
   predicted_ranks_l <- as.numeric(predicted_ranks_l)
   cluster_probs_l <- predict(mymodel, ndata.kd, type = 'probs') #cluster_probs
   
   for (i in seq(1, dim(cluster_probs_l)[1])){ 
      x <- as.numeric(cluster_probs_l[i,])
      max_1 <- max(x)
      max_1_pos <- which(x==max_1) 
      x[max_1_pos] <- 0
      max_2 <- max(x)
      if(abs(max_1-max_2) < PROB_DIFF_THRESHOLD) {
         predicted_ranks_l[i] <- 0
      }
      #break() 
   }
   return(predicted_ranks_l)
}


obtain_data_by_cluster <- function(data.sim, clusterCut){
   data_by_cluster <- list() 
   for(cluster_no in 1:length(unique(clusterCut))){
      idx <- paste('CL', cluster_no, sep = '')
      data_by_cluster[[idx]] <- data.sim[, clusterCut==cluster_no]
   }   
   return(data_by_cluster)
}


obtain_clusterCut <- function(ndata, num_clusters){ 
   # ndata: rows are observations, columns are variables
   # corMethod <- 'pearson' # 'spearman' # 
   # refCor <- cor((ndata), method = corMethod)
   # data.dist <- as.dist((1-refCor)/2)
   
   data.dist <- dist(t(ndata), method="euclidean")
   #data.dist <- dist(t(ndata), method="manhattan")
   #data.dist <- dist(t(ndata), method="minkowski")
   # data.dist <- dist(t(ndata), method="maximum")
   #data.dist <- dist(t(ndata), method="canberra")
   
   hc.rows <- hclust(data.dist, "ward.D2") #cluster along rows
   
   cutree_ranks = cutree(hc.rows, k=num_clusters) 
   return(cutree_ranks)
} 


normalize_by_wt_mean_and_sd <- function(geneExpression, means, sds){
   geneExpression <- sweep(geneExpression, 1, means, FUN = "-")
   geneExpression <- sweep(geneExpression, 1, sds, FUN = "/")  
   return(geneExpression)
}

