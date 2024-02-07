
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


# obtain the maximim significance value for the circuit nodes 
# coming from Netact or MARINa
obtain.sigval.Netact_or_MARINa <- function(coreTFs.Method, circuit.nodes){
   # coreTFs.Method=coreTFs.RI
   # circuit.nodes=circuit.nodes
   circuit.nodesFromMethod <- coreTFs.Method[which(coreTFs.Method %in% circuit.nodes)]  
   if(length(circuit.nodesFromMethod)>0){
      # sigvals.Method <- c(min(as.numeric(labels(circuit.nodesFromMethod))), 
      #                     max(as.numeric(labels(circuit.nodesFromMethod)))) 
      maxSigval <- max(as.numeric(labels(circuit.nodesFromMethod))) 
   }else{
      maxSigval <- -100
   }
   return(maxSigval) 
}

obtain.sigval.RI <- function(coreTFs.Method, circuit.nodes){
   # coreTFs.Method=coreTFs.RI
   # circuit.nodes=circuit.nodes
   circuit.nodesFromMethod <- coreTFs.Method[which(coreTFs.Method %in% circuit.nodes)]  
   if(length(circuit.nodesFromMethod)>0){
      # sigvals.Method <- c(min(as.numeric(labels(circuit.nodesFromMethod))), 
      #                     max(as.numeric(labels(circuit.nodesFromMethod)))) 
      minSigval <- min(as.numeric(labels(circuit.nodesFromMethod))) 
   }else{
      minSigval <- -100
   }
   return(minSigval) 
}





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


# Calculate average distance (sracipeHeatmapSimilarity) using boot strap sampling
cal.avgDist.bootStrap <- function(data.REF, data.sim, clusterCut.REF, BOOTSTRAP.ITER){
   AvgDist <- vector(mode = 'numeric', length = BOOTSTRAP.ITER)
   for(b in 1:BOOTSTRAP.ITER){
      print(b) 
      # calculate similarity between activitites and racipe simulation data:
      hS <- sracipeHeatmapSimilarity(dataReference = data.REF,
                                     #dataSimulation = data.sim, 
                                     dataSimulation = data.sim[, sample(1:dim(data.sim)[2], replace = TRUE)], # bootstrap sampling
                                     returnData = F, 
                                     #nClusters = NUM_CLUSTERS, 
                                     clusterCut = clusterCut.REF)
      AvgDist[b] <- hS$AvgDist 
      #break()
   } 
   return(AvgDist)
}

# Three way sorting
#-----------------
# Sort features of the simulated circuits:
# (1) Sort data by Accuracy (decreasing order) and add index to the data  
# (2) Sort the sorted data in step 1 by AvgDist (increading order) and add another index 
# (3) Add two indices to obtain the conjugate indix 
# (4) Sort the data from step 3 in increasing order (decreasing order) 
# (5) return sorted data from step 5
sortByTwoIndices <- function(circuit_metrics.sim){
   # sort by Accuracy:
   data.sortedByAccuracy <- circuit_metrics.sim[order(circuit_metrics.sim$Accuracy, decreasing = T), ]
   data.sortedByAccuracy$idxAccuracy <- 1:dim(data.sortedByAccuracy)[1]
   
   # sort by AvgDist:
   data.sortedByAvgDist <- data.sortedByAccuracy[order(data.sortedByAccuracy$AvgDist, decreasing = F), ]
   data.sortedByAvgDist$idxAvgDist <- 1:dim(data.sortedByAvgDist)[1]
   
   # add both indices:
   data.sortedByAvgDist$idxBoth <- data.sortedByAvgDist$idxAccuracy+data.sortedByAvgDist$idxAvgDist
   
   # sort by conjugate index:
   data.sortedByBoth <- data.sortedByAvgDist[order(data.sortedByAvgDist$idxBoth, decreasing = F),]   
   return(data.sortedByBoth)
}


# sort by accuracy and flexibility
sortByTwoIndices.acc_and_flex <- function(circuit_metrics.sim){
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




# sort by avgDistance and flexibility
sortByTwoIndices.avgDist_and_flex <- function(circuit_metrics.sim){
   # sort by Accuracy:
   # data.sortedByAccuracy <- circuit_metrics.sim[order(circuit_metrics.sim$Accuracy, decreasing = T), ]
   # data.sortedByAccuracy$idxAccuracy <- 1:dim(data.sortedByAccuracy)[1]
   
   # # sort by AvgDist:
   data.sortedByAvgDist <- circuit_metrics.sim[order(circuit_metrics.sim$AvgDist, decreasing = F), ]
   data.sortedByAvgDist$idxAvgDist <- 1:dim(data.sortedByAvgDist)[1]
   
   # sort by flexibility:
   data.sortedByFlexibility <- data.sortedByAvgDist[order(data.sortedByAvgDist$flexibility, decreasing = T), ]
   data.sortedByFlexibility$idxFlexibility <- 1:dim(data.sortedByFlexibility)[1] # adding index
   
   # add both indices: 
   data.sortedByFlexibility$idxBoth <- data.sortedByFlexibility$idxAvgDist+data.sortedByFlexibility$idxFlexibility
   
   # sort by conjugate index:
   data.sortedByBoth <- data.sortedByFlexibility[order(data.sortedByFlexibility$idxBoth, decreasing = F),]   
   return(data.sortedByBoth)
}


# Four way sorting
#-----------------
# Sort features of the simulated circuits:
# (1) Sort data by Accuracy (decreasing order) and add index to the data  
# (2) Sort the sorted data in step 1 by AvgDist (increading order) and add another index 
# (3) Add two indices to obtain the conjugate indix 
# (4) Sort the data from step 3 in increasing order (decreasing order) 
# (5) return sorted data from step 5
sortByThreeIndices <- function(circuit_metrics.sim){
   # sort by Accuracy:
   data.sortedByAccuracy <- circuit_metrics.sim[order(circuit_metrics.sim$Accuracy, decreasing = T), ]
   data.sortedByAccuracy$idxAccuracy <- 1:dim(data.sortedByAccuracy)[1] # adding index
   
   # sort by AvgDist:
   data.sortedByAvgDist <- data.sortedByAccuracy[order(data.sortedByAccuracy$AvgDist, decreasing = F), ]
   data.sortedByAvgDist$idxAvgDist <- 1:dim(data.sortedByAvgDist)[1] 
   
   # add both indices:
   data.sortedByAvgDist$idxBoth <- data.sortedByAvgDist$idxAccuracy+data.sortedByAvgDist$idxAvgDist
   
   # sort by flexibility:
   data.sortedByFlexibility <- data.sortedByAvgDist[order(data.sortedByAvgDist$flexibility, decreasing = T), ]
   data.sortedByFlexibility$idxFlexibility <- 1:dim(data.sortedByFlexibility)[1] # adding index
   
   
   # add three indices:
   data.sortedByFlexibility$idxTrio <- data.sortedByFlexibility$idxAccuracy+data.sortedByFlexibility$idxAvgDist + data.sortedByFlexibility$idxFlexibility
   
   
   # sort by conjugate index:
   data.sortedByTrio <- data.sortedByFlexibility[order(data.sortedByFlexibility$idxTrio, decreasing = F),]   
   return(data.sortedByTrio)
}



# Calculates specific features of a list of circuits 
# indexed by a conjugate key "no_top_tfs"-"feature_ratio_cutoff"-"absCor"
cal.circuit_metrics.globally <- function(networkList){ 
   library(igraph)
   # create a matrix to save features of each circuit: 
   circuit.features <- c('FeatureRatio', 'TopTFs', 'AbsCor', 
                         'Nodes', 'Interactions', 'PosInt', 
                         'Connected', 'Transitivity', 'MeanDistance')
   
   networkMetrics <- as.data.frame(matrix(nrow = length(circuits), ncol = length(circuit.features)))
   colnames(networkMetrics) <-  circuit.features   
   rownames(networkMetrics) <- names(circuits) 
   
   # loop through each circuit and calculate the features of the circuit:
   for(circuit_idx in names(networkList)){
      feature_ratio_cutoff  <- strsplit(circuit_idx, '-', 2)[[1]][1] 
      no_top_tfs <- strsplit(circuit_idx, '-', 2)[[1]][2] 
      absCor <- strsplit(circuit_idx, '-', 2)[[1]][3] 
      
      circuit <- networkList[[circuit_idx]]
      
      g <- graph_from_data_frame(circuit, 
                                 directed = TRUE, 
                                 vertices = NULL) 
      stat.vector <- c(feature_ratio_cutoff = as.numeric(feature_ratio_cutoff),
                       no_top_tfs = as.numeric(no_top_tfs),  
                       absCor = as.numeric(absCor), 
                       Nodes = length(union(circuit$Source,circuit$Target)),
                       Interactions = length(circuit$Source),
                       PosInt = length(which(circuit$Type ==1)),
                       Connected = igraph::is_connected(g), 
                       Transitivity = igraph::transitivity(g),
                       MeanDistance = igraph::mean_distance(g, directed = TRUE, 
                                                            unconnected = FALSE) 
      )
      networkMetrics[circuit_idx, ] <- stat.vector
   } 
   return(networkMetrics)
}


# find the duplicated topologies in list of circuits
find.tpo_dup_status <- function(circuit.list){
   network.similarity <- cal.network_metrics.pre_sim(networkList = circuit.list) 
   network.similarity.tmp <- network.similarity
   network.similarity.sorted <- network.similarity[with(network.similarity, order(Nodes, Interactions, PosInt)), ] 
   network.similarity <- network.similarity.sorted
   
   # Find duplicate topologies 
   # screening by metric node_count 
   node_count.similarity <- duplicated(network.similarity$Nodes)
   sum(node_count.similarity)   
   node_count.similarity   
   
   # screening by metric interaction_count.similarity
   interaction_count.similarity <- duplicated(network.similarity$Interactions) 
   interaction_count.similarity
   sum(interaction_count.similarity)
   
   # screening by metric posInt.similarity
   posInt.similarity <- duplicated(network.similarity$PosInt)  
   sum(posInt.similarity)
   posInt.similarity
   
   # find duplicates
   node_interaction_posInt.similarity <- node_count.similarity & interaction_count.similarity & posInt.similarity
   sum(node_interaction_posInt.similarity)
   
   network.similarity$DUP.STATUS <- node_interaction_posInt.similarity   
   return(network.similarity)
}


# Calculate adjacency matrix from a given topology
#  
cal.adj_matrix_from_tpo <- function(net_tpo){
   nodes <- sort(union(net_tpo$Source, net_tpo$Target))
   adj_matrix <- matrix(nrow = length(nodes), ncol = length(nodes)) 
   dim(adj_matrix)
   colnames(adj_matrix) <- nodes
   rownames(adj_matrix) <- nodes
   adj_matrix[] <- 0L
   
   for(row.count in 1:dim(net_tpo)[1]){ 
      node.source <- as.character(net_tpo$Source[row.count])
      node.target <- as.character(net_tpo$Target[row.count])
      adj_matrix[node.source, node.target] <- as.numeric(net_tpo$Type[row.count])
   }   
   return(adj_matrix)
} 


# Infer circuits for a set of core TFs based on a given 
# TF-target DB
infer.circuits.by.targetDB <- function(coreTFs, 
                                       targetDB, 
                                       eset.brain_array, 
                                       de.results, 
                                       int.strengths, 
                                       SUBNETWORK.SIZE.TSH){ 
   ## parameters:
   # coreTFs = coreTFs
   # targetDB = targetDB.list[[fr]]  
   # eset.brain_array = eset.brain_array
   # de.results = de.results
   # int.strengths <- INTERACTION.STRENGTHS
   
   ## databases and data structures of this function: 
   # coreTFs
   # targetDB 
   # eset.brain_array 
   # de.results  
   # functions: 
   # createNetworkMiValue
   # retain_uniq_networks 
   # cal.network_metrics.pre_sim 
   # select.largest.connected.subgraph 
   
   # Construct initial network 
   # (1) Construct initial network from TF-target DB
   #------------------------------------------------
   network <- data.frame(Source = character(), Target = character())
   network <- data.frame(Source = rep(names(lapply(targetDB, function(x)names(x))),
                                      times = unlist(lapply(targetDB, function(x)length(x)))),
                         Target = unlist(targetDB))
   
   
   # (2) Retain subnetwork restricted to only core TFs
   #--------------------------------------------------
   # retain TFs in TF-target DB
   #---------------------------
   length(coreTFs)
   coreTFs <- intersect(coreTFs, names(targetDB)) 
   length(coreTFs)
   
   # Retain subnetwork restricted to only core TFs
   #----------------------------------------------
   network <- network[intersect(which(network$Source %in% coreTFs),
                                which(network$Target %in% coreTFs)),] 
   
   # (3) Retain subnetwork restricted to activities 
   # Calculate TF activities  
   #------------------------
   a = TF_Activity(tfs = coreTFs,
                   GSDB = targetDB,  
                   eset = eset.brain_array,
                   DErslt = de.results  #DErslt=de.results$Overall
   ) 
   TFactivities <- a$all_activities #a.activities
   
   
   # Retain subnetwork restricted to activities 
   #---------------------------------------------------------
   network <- network[intersect(which(network$Source %in% rownames(TFactivities)),
                                which(network$Target %in% rownames(TFactivities))),]
   
   
   # (4) Add interactions (+ve or -ve type) to the network
   #------------------------------------------------------
   # find network genes
   networkGenes <- union(network$Source, network$Target)
   length(networkGenes)
   
   # subset expression data related to network genes
   TFactivities <- t(TFactivities[networkGenes,])
   dim(TFactivities)
   
   # Infer and assign interaction types
   #------------------------------------ 
   # calculate gene expression correlations
   actCor <- cor(TFactivities, method = "s") 
   dim(actCor)
   #head(rownames(actCor))
   #head(actCor, 2)
   
   # calculate correlation between Source and Target
   int.type <- integer(length = length(network$Source))
   for(i in seq_along(network$Source)){
      int.type[i] <- actCor[as.character(network$Source[i]), as.character(network$Target[i])]
   }
   
   # correlation
   network$Cor <- int.type
   
   # Deterimine interaction type from the Source:Target correlations
   # Assign interaction sign based on correlation
   # Convert correlations to interaction type 
   # >0: 1 (activation)
   # <0: 2 (inhibition)
   int.type[int.type>0] <- 1
   int.type[int.type<=0] <- 2
   network$Type <- int.type
   
   # assign absolute value of correlations to MI
   network$Mi <- abs(network$Cor)
   
   
   # Infer CANDIDATE networks from initial network 
   #----------------------------------------------
   print("Infer CANDIDATE networks")
   networkList  <- list() 
   for(int.strength in int.strengths){
      networkList[[as.character(int.strength)]] <- createNetworkMiValue(tmpNetwork = network, 
                                                                        coreTf = networkGenes, 
                                                                        posMi = int.strength, 
                                                                        negMi = int.strength, 
                                                                        name="net", 
                                                                        simulate = F)
   } 
   
   #names(networkList)
   #length(names(networkList))
   
   
   # Remove duplicated networks based on number of nodes and number of interactions 
   # (apply only once)
   #--------------------------------------------------------------------------------
   #length(networkList)
   networkList <- retain_uniq_networks(networkList)
   #length(networkList)
   
   
   # Refine sampled networks: 
   # Retain largest connected component having
   # nodes more than SUBNETWORK.SIZE.TSH  
   ##----------------------------------------- 
   networkList.refined <- list() 
   for(int.strength in names(networkList)){
      print(int.strength) 
      network.cur <- networkList[[int.strength]]  
      network.comp <- select.largest.connected.subgraph(network.cur=network.cur, 
                                                        network = network, 
                                                        subnetwork.size.TSH = SUBNETWORK.SIZE.TSH)
      if(!is.null(network.comp))
         networkList.refined[[int.strength]] <- network.comp
   }
   
   names(networkList.refined)
   length(names(networkList.refined))
   
   return(networkList.refined)
}


# Select the TOP n number of TFs from the combined TFs 
# from each of three methods: NetAct, MARINa, and RI
sel.topTFs.from_each_method <- function(TFs.combined, 
                                        NO_TOP_TFS_FROM_EACH_METHOD = 30){
   # TFs from netact
   #----------------
   tfs.netact <- TFs.combined$NetAct
   length(tfs.netact$tf)  
   
   # order the data frame by increasing value of z and q values:  
   tfs.netact <- tfs.netact[order(tfs.netact$zq, decreasing = FALSE), ]
   
   
   # TFs from MARINa
   #----------------
   tfs.MARINa <- TFs.combined$MARINa
   # order the data frame by increasing value of FDR 
   tfs.MARINa  <- tfs.MARINa[order(tfs.MARINa$FDR, decreasing = FALSE), ] 
   
   
   # TFs deom RI method
   #---------------------
   tfs.RI <- TFs.combined$RI
   colnames(tfs.RI)  
   # order the data frame by decreasing value of aggregate error:
   tfs.RI <- tfs.RI[order(tfs.RI$aggr.error, decreasing = TRUE), ]
   
   
   # Select top TFs from netact method: 
   if(length(tfs.netact$tf)<=NO_TOP_TFS_FROM_EACH_METHOD){
      top.tfs.netact <- tfs.netact$tf
   } else {
      top.tfs.netact <- tfs.netact$tf[1:NO_TOP_TFS_FROM_EACH_METHOD] 
   }
   
   # Select top TFs from MARINa method: 
   if(length(tfs.MARINa$tf)<=NO_TOP_TFS_FROM_EACH_METHOD){
      top.tfs.MARINa <- tfs.MARINa$tf
   } else {
      top.tfs.MARINa <- tfs.MARINa$tf[1:NO_TOP_TFS_FROM_EACH_METHOD] 
   }
   
   # Select top TFs from RI method: 
   if(length(tfs.RI$tf)<=NO_TOP_TFS_FROM_EACH_METHOD){
      top.tfs.RI <- tfs.RI$tf
   } else {
      top.tfs.RI <- tfs.RI$tf[1:NO_TOP_TFS_FROM_EACH_METHOD] 
   }
   
   # Combine TOP TFs from all THREE METHODS
   tfs.COM <- union(top.tfs.netact, 
                    union(top.tfs.MARINa, top.tfs.RI)
   )
   
   length(tfs.COM) # ==> 76
   length(unique(tfs.COM))    
   
   return(tfs.COM)
}



select.largest.connected.subgraph <- function(network.cur,
                                              network,
                                              subnetwork.size.TSH = 0.80){
   #node.total <- union(network$Source, network$Target) 
   node.total <- union(network.cur$Source, network.cur$Target)
   length(node.total)
   g <- igraph::graph_from_data_frame(network.cur, 
                                      directed = TRUE, 
                                      vertices = NULL)    
   clu <- igraph::components(g)
   idx <- which(clu$csize==max(clu$csize))
   
   g.groups <- igraph::groups(clu) 
   
   # if there are more than one subnetworks of the largest size, select the first one:
   nodes.largest_component <- g.groups[[idx[1]]] 
   
   #prop.nodes <- length(nodes.largest_component)/length(node.total)
   
   if(length(nodes.largest_component)/length(node.total) <= subnetwork.size.TSH){
      return() #return(NULL)
   }
   
   status.SOURCE <- network$Source %in% nodes.largest_component
   status.TARGET <- network$Target %in% nodes.largest_component
   status.BOTH <- status.SOURCE & status.TARGET
   network.sele <- network[status.BOTH,] 
   return(network.sele)
}



cal.network_metrics.pre_sim <- function(networkList){
   networkMetrics <- NULL
   for(mi in names(networkList)){
      #print(mi) 
      circuit <- networkList[[mi]]
      
      posMi <- as.numeric(mi)
      negMi <- as.numeric(mi)
      
      #print(dim(circuit))
      g <- graph_from_data_frame(circuit, 
                                 directed = TRUE, 
                                 vertices = NULL)
      networkMetricsTmp <- data.frame(MI=as.numeric(mi), 
                                      Nodes= length(union(circuit$Source,circuit$Target)),
                                      Interactions=length(circuit$Source),
                                      PosInt = length(which(circuit$Type ==1)),
                                      Connected = igraph::is_connected(g), 
                                      Transitivity = igraph::transitivity(g),
                                      MeanDistance = igraph::mean_distance(g, directed = TRUE, 
                                                                           unconnected = FALSE), 
                                      MiPos = posMi, 
                                      MiNeg = negMi)  
      networkMetrics <- rbind(networkMetrics, networkMetricsTmp)
   }   
   return(networkMetrics)
}

# Remove duplicated networks based on number of nodes and number of interactions
#------------------------------------------------------------------------------
retain_uniq_networks <- function(networkList){
   # find unique networks
   network.similarity <- cal.network_metrics.pre_sim(networkList = networkList) 
   
   # sort the networks by three network properties:
   network.similarity.sorted <- network.similarity[with(network.similarity, order(Nodes, Interactions, PosInt)), ] 
   network.similarity <- network.similarity.sorted
   
   node_count.similarity <- duplicated(network.similarity$Nodes)
   sum(node_count.similarity)   
   node_count.similarity   
   interaction_count.similarity <- duplicated(network.similarity$Interactions)
   interaction_count.similarity
   sum(interaction_count.similarity)
   
   posInt.similarity <- duplicated(network.similarity$PosInt) 
   
   # find duplicated networks:
   node_interaction.similarity <- node_count.similarity & interaction_count.similarity & posInt.similarity
   sum(node_interaction.similarity)

   # retain the unique networks 
   MIs.retained <- network.similarity$MI[!node_interaction.similarity]
   MIs.retained <- unlist(lapply(MIs.retained, function (x) as.character(x)))
   networkList <-  networkList[MIs.retained]   
   return(networkList)
}


cal.network_metrics.post_sim <- function(simulatedNetworks, expData.REFERENCE){
   #library(sRACIPE) 
   
   networkMetrics <- NULL
   nClusters <- 2
   for(mi in names(simulatedNetworks)){    
      posMi <- as.numeric(mi)
      negMi <- as.numeric(mi)
      
      racipe <- simulatedNetworks[[mi]]
      
      #racipe <- sracipeNormalize(racipe)
      expData <- assay(racipe,1)
      
      tmp1 <- expData.REFERENCE 
      tmp2 <-  expData
      nClusters <- 2
      
      hS <- sracipeHeatmapSimilarity(dataReference = tmp1,
                                     dataSimulation = tmp2,
                                     nClusters = nClusters) 
      
      ClusterA <- hS$simulated.cluster.freq[2]
      ClusterB <- hS$simulated.cluster.freq[3]
      
      circuit <- as.data.frame(sracipeCircuit(racipe)) 
      g <- graph_from_data_frame(circuit, directed = TRUE, vertices = NULL)
      expDataTmp <- as.matrix(expData)
      storage.mode(expDataTmp) <- "logical"
      networkMetricsTmp <- data.frame(MI=as.numeric(mi), 
                                      Nodes= length(union(circuit$Source,circuit$Target)),
                                      Interactions=length(circuit$Source),
                                      PosInt = length(which(circuit$Type ==1)),
                                      ClusterA = ClusterA, ClusterB = ClusterB,
                                      Connected = igraph::is_connected(g), 
                                      ModelCount = sracipeConfig(racipe)$simParams["numModels"], 
                                      Transitivity = igraph::transitivity(g),
                                      MeanDistance = igraph::mean_distance(g, directed = TRUE, 
                                                                           unconnected = FALSE),
                                      MiPos = posMi, MiNeg = negMi)
      networkMetricsTmp$Accuracy <- (networkMetricsTmp$ClusterA + networkMetricsTmp$ClusterB)
      networkMetricsTmp$Accuracy2 <- hS$KL
      networkMetricsTmp$ClusterA2 <- hS$cluster.similarity[2]
      networkMetricsTmp$ClusterB2 <- hS$cluster.similarity[3]
      networkMetrics <- rbind(networkMetrics, networkMetricsTmp) 
   }   
   return(networkMetrics)
}





# Function to generate the networks  
#==================================== 
#' @param tmpNetwork A network containing the TF, target and the metric.
#' The metric can be mutual information. It should also have a metric to 
#' identify the sign of the interaction, for example, correlation.
#' In this case the network can have same interaction with different
#' mutual information and/or sign metric as the interaction can come from
#' different sources.
#' @param coreTf List of factors whose nearest neighbors are to be included in 
#' the network. Use all tfs as coreTfs if no core tf is to be used.
#' @param posMi Cutoff for positive mutual information
#' @param negMi Cutoff for negative mutual information
#' @param name generic name for dataset
#' @param simulate whether to simulate the networks using sRACIPE
#' @param allowAllSameSign wether to simulate networks with interactions of
#' one type only
#' @param minCoreTFs minimum number of coreTFs that should be in the network
#' for it to be included in simulations
#' @param minPosInt minimum number of excitatory interactions that should be in the network
#' for it to be included in simulations.
#' @param minNegInt minimum number of inhibitory interactions that should be in the network
#' for it to be included in simulations.
#'
#' 
createNetworkMiValue <- function(tmpNetwork=network, 
                                 coreTf = coreTf, 
                                 posMi = posMi,
                                 negMi = negMi, name = "A549", 
                                 simulate = FALSE, allowAllSameSign = FALSE,
                                 minCoreTFs = 0, minPosInt=0,minNegInt=0){
   # tmpNetwork <- tmpNetwork[union(which(tmpNetwork$Source %in% coreTf),
   #                                which(tmpNetwork$Source %in% coreTf)),]
   # sort the network based on mutual information
   tmpNetwork <- tmpNetwork[order(tmpNetwork$Mi, decreasing = TRUE),]
   # remove any duplicate interactions
   tmpNetwork <- tmpNetwork[!duplicated(tmpNetwork[c("Source","Target","Type")]),]
   # Select the positive network
   posNetwork <- tmpNetwork[tmpNetwork$Cor >0,]
   posNetwork <- posNetwork[which(posNetwork$Mi > posMi),]
   # Select the negative network
   negNetwork <- tmpNetwork[tmpNetwork$Cor < 0,]
   negNetwork <- negNetwork[which(negNetwork$Mi > negMi),]
   
   # Exit if there are either no positive or negative interaction.
   # Can be removed if all positive or all negative networks are ok
   if(!allowAllSameSign){
      if(dim(posNetwork)[1] == 0 | dim(negNetwork)[1] == 0) return()
   }
   # Combine the network
   tmpNetwork <- rbind(posNetwork,negNetwork)
   print("check if there are any duplicate interactions with same or conflicting sign")
   print(sum(duplicated(tmpNetwork[,c("Source","Target")])))
   print(" check if there are any duplicate interactions with same sign")
   print(sum(duplicated(tmpNetwork[,c("Source","Target","Type")])))
   # Remove duplicate interactions
   tmpNetwork <- tmpNetwork[!duplicated(tmpNetwork[,c("Source","Target","Type")]),]
   print("Check if any conflicting interaction is there")
   print(tmpNetwork[duplicated(tmpNetwork[,c("Source","Target")]),])
   print(tmpNetwork[duplicated(tmpNetwork[,c("Source","Target")], fromLast = TRUE),])
   # Remove conflicting interactions
   delRows <- anyDuplicated(tmpNetwork[,c("Source","Target")])
   delRows <- c(delRows, anyDuplicated(tmpNetwork[,c("Source","Target")], fromLast = TRUE))
   delRows <- delRows[delRows>0]
   if(length(delRows>0)){ tmpNetwork <- tmpNetwork[-delRows,]}
   
   # Remove signaling interactions
   # Uncomment for iterative removal
   # interactionRemoved = length(tmpNetwork$Source)
   # while(interactionRemoved>0)
   {
      tmpVar = length(tmpNetwork$Source)
      tmpNetwork <- tmpNetwork[(which(tmpNetwork$Source %in% tmpNetwork$Target )),]
      # if targets only nodes are also to be removed.
      # tmpNetwork <- tmpNetwork[which(tmpNetwork$Target %in% tmpNetwork$Source),]
      # interactionRemoved = tmpVar - length(tmpNetwork$Source)
   }
   
   tmpNetwork <- tmpNetwork[!duplicated(tmpNetwork),]
   #tmpNetwork <- removeNodePairs(tmpNetwork) 
   print("Number of interactions")
   print(length(tmpNetwork$Source))
   networkTfs <- union(tmpNetwork$Source,tmpNetwork$Target)
   #print("Core Tfs included in the network")
   #print(networkTfs[which(networkTfs %in% coreTfGenes)])
   
   #if(length(networkTfs) < minCoreTFs) simulate = FALSE
   
   if(length(which(tmpNetwork$Type ==2)) <minNegInt) simulate = FALSE
   if(length(which(tmpNetwork$Type ==1)) <minNegInt) simulate = FALSE
   require(igraph)
   circuit <-  tmpNetwork[,c("Source", "Target", "Type")]
   g <- graph_from_data_frame(circuit, directed = TRUE, vertices = NULL)
   if(!igraph::is_connected(g)) simulate = FALSE 
   
   
   if(simulate){
      require(sRACIPE)
      tmp <- sracipeSimulate(circuit = tmpNetwork[,c("Source", "Target", "Type")], plotToFile = TRUE,
                             numModels = 2000, plots = FALSE, stepper = "RK4", integrateStepSize = 0.05)
      annotation(tmp) <- paste0(name,"_",as.character(posMi), "_",as.character(negMi),"Interactions")
      # tmp <- sracipePlotData(tmp, plotToFile = T)
      saveRDS(tmp, file = paste0("simAll/",name,"_",as.character(posMi), "_",as.character(negMi),".rds"))
   }
   return(tmpNetwork) 
}


# Helper functions 
#-----------------
removeNodePairs <- function(tmpNetwork){
   srcGenes <- unique(tmpNetwork$Source)
   for(i in seq_along(srcGenes)){
      src <- tmpNetwork$Source[which(tmpNetwork$Source == srcGenes[i])]
      tgtGenes <- tmpNetwork$Target[which(tmpNetwork$Target == srcGenes[i])]
      if((length(tgtGenes)==1) & (length(src)==1) ){
         if((tgtGenes==src)){
            print(tgtGenes)
            tmpNetwork <- tmpNetwork[-(which(tmpNetwork$Source == tgtGenes)),]
            tmpNetwork <- tmpNetwork[-(which(tmpNetwork$Target == tgtGenes)),]
         }
      }
   }
   return(tmpNetwork)
}




show_network.dynamically <- function(top.df) {
   # top.df - column 1: SOURCE, column 2: TARGET, column 3: TYPE
   library(visNetwork)
   #n.list <- as.character(unique(union(top.df$SOURCE, top.df$TARGET)))
   n.list <- as.character(unique(union(top.df[,1], top.df[,2])))
   nodes.df <- data.frame(id = n.list,
                          label =  n.list,
                          font.size =30,
                          shape='circle')

   edge_col.df <- data.frame(c(1,2),c("blue","darkred"))
   colnames(edge_col.df) <- c("interaction", "color")

   #color=edge_col.df$color[top.df$TYPE]
   color=edge_col.df$color[top.df[,3]]

   #edges.df <- data.frame(from = as.character(top.df$SOURCE),
   #                       to = as.character(top.df$TARGET),
   #                       color=color)

   edges.df <- data.frame(from = as.character(top.df[,1]),
                          to = as.character(top.df[,2]),
                          color=color)


   visNetwork(nodes=nodes.df, edges = edges.df,
              height = "500px", width = "500%") %>%
      visEdges(arrows = "to") %>%
      visOptions(manipulation = TRUE) %>%
      visLayout(randomSeed = 123) %>%
      visPhysics(solver = "forceAtlas2Based")
}


#---------------------------------------------------------------------#
#                       map_probes_to_genes(probes)
#        maps a vector of probes to a data frame of mapped genes
#---------------------------------------------------------------------#
"
Input: a vector of probe names. 
Output:a dataframe with three columns: 
one: ENTREZID
second: SYMBOL 
third: GENENAME 
each row of the output dataframe corresponds to a probe 
in the input vector of the probe names in the same 
order. 
If no mapping is found, NA is used in the mapped name.
"
map_probes_to_genes <- function(probes_p){
   library(biomaRt)
   library(org.Hs.eg.db)
   
   probes.sel <- sapply(probes_p, 
                        function(dat) {strsplit(as.character(dat),
                                                ".",fixed=T)[[1]][1]})
   probes.1 <- as.character(probes.sel) 
   symbols <- mget(probes.1,org.Hs.egREFSEQ2EG,ifnotfound=NA)
   entre_ids <- as.character(symbols)
   mapped_genes <- biomaRt::select(org.Hs.eg.db, entre_ids, 
                                   c("SYMBOL","ENTREZID", "GENENAME"))
   return(mapped_genes)
}


