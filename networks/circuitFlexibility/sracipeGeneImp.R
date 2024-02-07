#' @export
#' @import philentropy Hmisc
#' @title  A method to evaluate the effects of gene perturbation without clustering
#' @description Computes the distance between distribution of principle components
#' of the original data and subsets of data obtained by selecting the models
#' for which production rate of genes are in bottom kdPercent for knockdown and 
#' more than oePercent for overexpression. Distance of each PC is weighted by
#' the variance explained by that PC.
#' @param racipe RacipeSE object
#' @param kdPercent numeric Default 10. Below what percent should be the 
#' production rate of a gene for selection of models for knockdown.
#' @param oePercent numeric Default 90. Above what percent should be the 
#' production rate of a gene for selection of models for overexpression.
#' @param nBins Number of bins used for binning the data.
#' @param distMethod Distance method to be used for calculating the distance. 
#' List of implemented probability distance measures can be obtained using
#' philentropy::getDistMethods()
#' @param ... Additional arguments, for use in specific methods.
#' @examples
#' 
#' data("demoCircuit")
#' rSet <- sRACIPE::sracipeSimulate(circuit = demoCircuit, numModels = 100)
#' sracipeGeneImp(rSet)
#' @return data.frame A dataframe containing the distance for each gene 
#' perturbation (knockdown and overexpression). Using colMeans() on the result
#' provides an estimate of the rigidity of the network.
#'
sracipeGeneImp <- function(racipe = racipe, kdPercent=10, oePercent=90, 
                           nBins=12, distMethod =  "bhattacharyya",...
){
  require("philentropy")
  require("Hmisc")
  if(missing(nBins)) nBins = round(1+log2(metadata(racipe)$config$simParams["numModels"])) # Sturges' rule 
  if(!metadata(racipe)$normalized) racipe <- sracipeNormalize(racipe)
  expMat <- assay(racipe,1)
  pca <- summary(prcomp(t(expMat)))
  params <- sracipeParams(racipe)
  geneNames <- names(racipe)
  res <- data.frame(row.names = geneNames)
  kdSelf <- vector(mode = "numeric", length = length(geneNames))
  oeSelf <- vector(mode = "numeric", length = length(geneNames))
  for(i in seq_along(geneNames)){
    selParams <- params[paste0("G_",geneNames[i])]
    selModels1 <- which(selParams <(min(selParams) + 
                                      0.01*kdPercent*(max(selParams)
                                                      -min(selParams))))
    selModels2 <- which(selParams >(min(selParams) + 
                                      0.01*oePercent*(max(selParams) - min(selParams))))
    distPC <- 0.0
    for(pcCounter in seq_along(geneNames)){
      PCCount1 <- pca$x[,pcCounter]
      PCCount2 <- pca$x[selModels1,pcCounter]
      PCCount3 <- pca$x[selModels2,pcCounter]
      
      brkPoint <- Hmisc::cut2(PCCount1, g=nBins, onlycuts = T)
      PCCount1 <- table(Hmisc::cut2(PCCount1, cuts = brkPoint, minmax = F))/length(PCCount1)
      PCCount2 <- table(Hmisc::cut2(PCCount2, cuts = brkPoint, minmax = F))/length(PCCount2)
      PCCount3 <- table(Hmisc::cut2(PCCount3, cuts = brkPoint, minmax = F))/length(PCCount3)
      
      distVal1 <- philentropy::distance(rbind(PCCount1,PCCount2), method = distMethod, ...)
      distVal2 <- philentropy::distance(rbind(PCCount1,PCCount3), method = distMethod, ...)
      kdSelf[i] <- kdSelf[i] + pca$importance[2,pcCounter]*distVal1
      oeSelf[i] <- oeSelf[i] + pca$importance[2,pcCounter]*distVal2
    }
  }
  res$kdSelf <- kdSelf
  res$oeSelf <- oeSelf
  # return(colMeans(res))
  return(res)
}
