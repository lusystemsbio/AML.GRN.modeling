normalize_by_wt_mean_and_sd <- function(geneExpression, means, sds){
  geneExpression <- sweep(geneExpression, 1, means, FUN = "-")
  geneExpression <- sweep(geneExpression, 1, sds, FUN = "/")  
  return(geneExpression)
}

