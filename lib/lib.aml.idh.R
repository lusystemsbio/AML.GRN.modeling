#===========================================================#
#                    LIBRARY FOR FUNTION DEFINITIONS  
#                          lib.aml.idh.R
#                        aml.idh/phase.III/lib
#===========================================================#


#7 Combine Activity and Expression Heatmap
Combine_heatmap = function(new_activity, eset){ 
   library(Biobase)
   require(ComplexHeatmap)
   require(circlize)
   if (is(eset, "ExpressionSet")){
      data = exprs(eset)
   }else{data = eset}
   H1 = Heatmap(row_norm(new_activity), col = colorRamp2(c(-2, 0, 2), c("green3", "white", "red")),
                cluster_columns = F, cluster_rows = T, show_row_dend = FALSE, name = "Activity", column_title = "Activity")
   gs = rownames(new_activity)
   gc = colnames(new_activity)
   H2 = Heatmap(row_norm(data[gs, gc]), col = colorRamp2(c(-2, 0, 2), c("green3", "white", "red")),
                cluster_columns = F, cluster_rows = T, show_row_dend = FALSE,name = "Expression", column_title = "Expression")
   #H1 + H2
   print(H1+H2)
}




#---------------------------------------------------------------------#
#                 Compute correlation distance 
#---------------------------------------------------------------------#
dist_cor <- function(x) as.dist(1-cor(t(x)))

cal_mean_sd <- function(data) {
   mean_sd_df <- as.data.frame(matrix(nrow = ncol(data),
                                      ncol = 2))
   NODE_LIST <- colnames(data)
   row.names(mean_sd_df) <- NODE_LIST
   names(mean_sd_df) <- c("MEAN", "SD")
   
   for (node in NODE_LIST) { 
      m <- mean(data[,node])
      sdt <- sd(data[,node])
      mean_sd_df[node,] <- c(m, sdt)
   }  
   return(mean_sd_df)
}

normalize_data <- function(data_l, mean_sd_df){
   ndata_l=data_l
   #for(i in 1:dim(data_l)[2]) 
   for(node in colnames(data_l))       
   {
      #print(node)
      ndata_l[,node] = (data_l[, node] - mean_sd_df[node,]$MEAN)/mean_sd_df[node,]$SD
   } 
   return(ndata_l)
}


#---------------------------------------------------------------#
# map.colnames_of_raw_data: maps column names of the raw AML data
#---------------------------------------------------------------#
map.colnames_of_raw_data <- function(fname_data, fname_mapping){
   #----------------------------------------------------------------#
   #                       Load data
   #----------------------------------------------------------------#
   cat("Loading expression data ...\n")
   # load raw data:
   rdata <- read.table(file = fname_data, header=TRUE, sep ='\t', row.names=1)
   
   cat("Loading mapping data ...\n")
   # load sample name mapping file: 
   map_data <- read.table(file = fname_mapping, header=TRUE, sep ='\t', row.names=NULL)
   
   label_vector <- map_data[,"label"]
   #t <- sapply(label_vector, function(k) paste0("D-", k))
   #colnames(rdata) <- sapply(label_vector, function(k) paste0("D-", k))
   
   # replace each '- in the column names with a '.'
   colnames(rdata) <- gsub('-', '.',colnames(rdata))
   
   class(map_data$file)
   
   cat("Mapping column names ...\n")
   
   # Create column names for CONTROL
   #================================== 
   
   # Extract file names and labels: 
   #-------------------------------
   ctrl.filenames <- as.character(map_data$file[map_data$type=='Control'])
   ctrl.labels <- as.character(map_data$label[map_data$type=='Control'])
   
   cat(ctrl.filenames)
   cat(ctrl.labels)
   
   # Create column names:
   #---------------------
   ctrl.labels <- sapply(ctrl.labels, function(k) paste0("CTRL.", k))
   cat(ctrl.labels)
   
   # replace - with . in the file names to match with colnames: 
   ctrl.filenames <- gsub('-', '.',ctrl.filenames)
   cat(ctrl.filenames)
   
   
   # Create column names for TREATMENTS
   #==================================== 
   
   # Extract file names and labels: 
   #-------------------------------
   trmt.filenames <- as.character(map_data$file[map_data$type=='AML'])
   trmt.labels <- as.character(map_data$label[map_data$type=='AML'])
   
   cat(trmt.filenames)
   cat("\n")
   cat(trmt.labels)
   cat("\n")
   # Create column names:
   #---------------------
   trmt.labels <- sapply(trmt.labels, function(k) paste0("AML.", k))
   cat(trmt.labels)
   
   # Replace column names with the newly created column names
   #========================================================= 
   library(data.table)
   
   # Replace column names for CONTROL
   #---------------------------------
   setnames(rdata, ctrl.filenames, ctrl.labels)
   
   # Replace column names for TREATMENT
   #-----------------------------------
   setnames(rdata, trmt.filenames, trmt.labels)
   
   cnames <- colnames(rdata)
   cat(cnames)
   return(rdata)
}


#---------------------------------------------------------------#
# bind.sample_names.3_groups: 
#   load and combine sample names
#   three groups: one CONTROL and two TREATMENT
#---------------------------------------------------------------#
bind.sample_names.3_groups <- function(samples.CTRL, 
                                       sample_ids.TREATMENT.1, 
                                       sample_ids.TREATMENT.2){
   
   # replace D with CTRL in sample names:
   #names(samples.CTRL) <- gsub('CONTROL_', 'CTRL.', colnames(samples.CTRL))
   sample_ids.CTRL <- as.character(unlist(samples.CTRL['SAMPLES',1:ncol(samples.CTRL)]))
   sample_ids.CTRL <- gsub('D', 'CTRL', sample_ids.CTRL)
   length(sample_ids.CTRL)
   cat(sample_ids.CTRL)
   
   # replace D- in the patient list with AML. : 
   sample_ids.TREATMENT.1 <- gsub('D-', 'AML.', sample_ids.TREATMENT.1)   
   length(sample_ids.TREATMENT.1)
   cat(sample_ids.TREATMENT.1)
   
   # replace D- in the patient list with AML. : 
   sample_ids.TREATMENT.2 <- gsub('D-', 'AML.', sample_ids.TREATMENT.2)   
   length(sample_ids.TREATMENT.2)
   cat(sample_ids.TREATMENT.2)
   
   # combine ids
   sample_ids.comb <- c(sample_ids.CTRL, 
                        sample_ids.TREATMENT.1, 
                        sample_ids.TREATMENT.2) 
   
   return(sample_ids.comb)
}


#---------------------------------------------------------------#
# create.colnames: 
#   load  
#   three groups: one CONTROL and two TREATMENT
#---------------------------------------------------------------#
create.colnames <- function(samples.CTRL, 
                            sample_ids.TREATMENT.1, 
                            sample_ids.TREATMENT.2, 
                            TREATMENT.1.NAME, 
                            TREATMENT.2.NAME) {
   
   # replace D with CTRL in sample names:
   #names(samples.CTRL) <- gsub('CONTROL_', 'CTRL.', colnames(samples.CTRL))
   sample_ids.CTRL <- as.character(unlist(samples.CTRL['SAMPLES',1:ncol(samples.CTRL)]))
   sample_ids.CTRL <- gsub('D', 'CTRL', sample_ids.CTRL)
   length(sample_ids.CTRL)
   cat(sample_ids.CTRL)
   
   # replace D- in the patient list with AML. : 
   sample_ids.TREATMENT.1 <- gsub('D-', paste0(TREATMENT.1.NAME, '.'), sample_ids.TREATMENT.1)   
   length(sample_ids.TREATMENT.1)
   cat(sample_ids.TREATMENT.1)
   
   # replace D- in the patient list with AML. : 
   sample_ids.TREATMENT.2 <- gsub('D-', paste0(TREATMENT.2.NAME, '.'), 
                                  sample_ids.TREATMENT.2)   
   length(sample_ids.TREATMENT.2)
   cat(sample_ids.TREATMENT.2)
   
   # combine ids
   colnames_for_data.sele <- c(sample_ids.CTRL, 
                        sample_ids.TREATMENT.1, 
                        sample_ids.TREATMENT.2)     
   
   return(colnames_for_data.sele)
}




#---------------------------------------------------------------#
# combine.sample_names.two_groups: 
#   load and combine sample names
#   two groups: one CONTROL and one TREATMENT
#---------------------------------------------------------------#
combine.sample_names.two_groups <- function(fname_samples.CTRL, 
                                              fname_samples.TREATMENT){
   
   # Load control ids:
   samples.CTRL <- read.csv(fname_samples.CTRL, row.names = 1)
   
   # replace D with CTRL in sample names:
   names(samples.CTRL) <- gsub('CONTROL_', 'CTRL.', colnames(samples.CTRL))
   sample_ids.CTRL <- as.character(unlist(samples.CTRL['SAMPLES',1:ncol(samples.CTRL)]))
   sample_ids.CTRL <- gsub('D', 'CTRL', sample_ids.CTRL)
   length(sample_ids.CTRL)
   cat(sample_ids.CTRL)
   
   # Load treatment ids: 
   samples.TREATMENT <- read.table(file = fname_samples.TREATMENT, 
                                   col.names = c("PATIENT.ID")) 
   
   # replace D- in the patient list with AML. : 
   sample_ids.TREATMENT <- gsub('D-', 'AML.', samples.TREATMENT$PATIENT.ID)   
   length(sample_ids.TREATMENT)
   cat(sample_ids.TREATMENT)
   
   # combine ids
   sample_ids.comb <- c(sample_ids.CTRL, sample_ids.TREATMENT)    
   return(sample_ids.comb)
}



get_new_colnames <- function(sample.names){ 
   #length(sample.names) # 39 
   # 11: CTRL 
   # 8: CL2B 
   # 10: CL3A
   # 4: CL4C
   # 6: CL5D
   samples.ctrl <- sample.names[1:11]
   samples.2b <- sample.names[12:19]
   samples.3a <- sample.names[20:29]
   samples.4c <- sample.names[30:33]
   samples.5d <- sample.names[34:39]
   samples.2b.new <- gsub('AML', 'CL2B', samples.2b)
   samples.3a.new <- gsub('AML', 'CL3A', samples.3a)
   samples.4c.new <- gsub('AML', 'CL4C', samples.4c) 
   samples.5d.new <- gsub('AML', 'CL5D', samples.5d) 
   
   sample.names.new <- c(samples.ctrl, samples.2b.new, samples.3a.new, 
                         samples.4c.new, samples.5d.new)
   
   return(sample.names.new)
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


#---------------------------------------------------------------------#
#                 create.ExpressionSet_instance (edata)
#  Creates an instance of ExpressionSet from the expression data
#---------------------------------------------------------------------#
create.ExpressionSet_instance <- function(edata){ 
   mapped.df <- map_probes_to_genes(row.names(edata)) # get probes to gene symbol mappings
   edata$Symbol <- mapped.df$SYMBOL 
   edata = edata[!is.na(edata[,"Symbol"]),]  # remove unannotated genes
   edata = aggregate(. ~ Symbol, edata, mean)  # mean expression for redudant probes
   rownames(edata) = edata$Symbol 
   edata$Symbol = NULL
   
   # convert edata to expression set type
   require(Biobase)
   eset <- new('ExpressionSet', exprs = as.matrix(edata))
   
   # construct phenoData
   #--------------------
   celltypes = c(rep("CTRL", 11), rep("IDH1", 10), rep("IDH2", 10))
   
   phenoData = new("AnnotatedDataFrame", data = data.frame(celltype = celltypes,  
                                                           samples = colnames(edata))) 
   rownames(phenoData) = colnames(edata)
   
   # Construct an instance of ExpressionSet 
   #---------------------------------------
   eset = ExpressionSet(assayData = as.matrix(edata), phenoData = phenoData) 
   
   return(eset)
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
#        cal.pairwise_mi <- function(data)
# Calculates pairwise MI of a set of random variables from 
# their observed values 
#---------------------------------------------------------------------#
# data: a matrix of observed values where
#     columns are random variables or gene names or nodes of a network
#     rows are conditions or observations 
#---------------------------------------------------------------------#
cal.pairwise_mi <- function(data, numBins1 = 8, numBins2 = 8){
   library(entropy) 
   
   nGenes <- dim(data)[2]

   # Create place holder for MI:
   mi_data <- matrix(0, nrow = nGenes, ncol = nGenes)
   geneNames <- colnames(data)
   rownames(mi_data) <- geneNames
   colnames(mi_data) <- geneNames
   
   # Calculate MI:
   for(i in 1:(nGenes-1)){
      for(j in (i+1):(nGenes)){
         temp <- entropy::discretize2d(data[,i], data[,j], 
                                       numBins1 = numBins1, 
                                       numBins2 = numBins2)
         mi_data[i, j] <- entropy::mi.shrink(temp, verbose = F)
         mi_data[j, i] <- mi_data[i, j]
      }
   } 
   diag(mi_data) <- 0
   return(mi_data)
}

#---------------------------------------------------------------------#
# build.network <- function(data, GSDB, MI_THRESHOLD)
# data: a matrix of observed values where
#     rows are conditions or observations 
#     columns are random variables or gene names or nodes of a network
# GSDB: transcription factor database 
# MI_THRESHOLD: MI threshold to be used to filter out TF-TF relations
#---------------------------------------------------------------------#
build.network <- function(data, mi_data, GSDB, MI_THRESHOLD){
   library(reshape2) 
   
   # Create TF-TF matrix with MI:
   miMat <- melt(mi_data)  
   
   # Compute pairwise correlations of TF:
   corMat <- cor(data) 
   dim(data)
   dim(corMat)

   # Select the source and target TFs whose pairwise MI is above threshold:    
   tf_source <- miMat[which(miMat$value > MI_THRESHOLD), 1] 
   tf_target <- miMat[which(miMat$value > MI_THRESHOLD), 2]

   # Construct SOURCE (TF) to TARGET (TF targets) mapping for the given gene set DB:
   netMat.gsdb <- allNet(GSDB) 
   
   # Select the entries from GSDB TF to target mappings where 
   #  both source and target belong to tf_source and tf_target, respectively 
   # and MI of each TF pair is above the MI threshold: 
   
   tf.network = data.frame(as.matrix(netMat.gsdb[netMat.gsdb$from %in% tf_source & 
                                                  netMat.gsdb$to %in% tf_target,]), 
                         stringsAsFactors = F)
   
   # Define the relation type of each TF-TF interaction pair:
   tf.network$relation <- NULL
   
   for (i in 1:length(tf.network$from)){
      source_name <- tf.network$from[i]
      target_name <- tf.network$to[i]
      
      # if the correlation is positive, the relation is excitatory, otherwise inhibitory:
      tf.network$relation[i] <- ifelse((corMat[source_name, target_name] > 0), 1, 2)
   }
   
   # set column names:
   colnames(tf.network) <- c("SOURCE", "TARGET", "TYPE")    
   
   return(tf.network)
}


#---------------------------------------------------------------------#
# order.TFs_with_Q0 (tfs.df)
# tfs.df: a data frame for TF and their ranked qvalues
#  return: Ordered TFs with zero Qvalue according to Z.score    
#---------------------------------------------------------------------#

order.TFs_with_Q0 <- function(tfs.df) {
   # Order TFs with Q-value Zero
   #------------------------------------------
   tfs.Q0.df <- tfs.df[tfs.df$qvals==0,]
   tfs.Qgt0.df <- tfs.df[tfs.df$qvals>0,]
   
   dim(tfs.Q0.df)
   dim(tfs.Qgt0.df)
   
   tfs.Q0.ordered.df <- tfs.Q0.df[with(tfs.Q0.df, order(-z)), ]
   
   #------------------------------------------
   # Add a new column to have both Z.score and Q values
   #------------------------------------------
   tfs.Q0.ordered.df$zq <- tfs.Q0.ordered.df$z
   tfs.Qgt0.df$zq <- tfs.Qgt0.df$qvals 
   
   #------------------------------------------
   # Put the two sub data frames together.
   #------------------------------------------
   tfs.ordered.df <- rbind(tfs.Q0.ordered.df, tfs.Qgt0.df)  
   
   return(tfs.ordered.df)
}




#---------------------------------------------------------------------#
# a NetAct function 
# Will not be needed after NetAct is updated 
#---------------------------------------------------------------------#
TF_Filter = function(actMat, GSDB, miTh = 1.4, nbins = 8, method = "spearman"){
   require(infotheo)
   #require(data.table)
   require(reshape2)
   NetMat = allNet(GSDB)
   # calculate the MI
   miMat = discretize(t(actMat), disc="equalfreq", nbins = 8)
   miMat = mutinformation(miMat, method = "shrink")
   diag(miMat) = 0
   actLinks = melt(miMat)
   corMat = cor(t(actMat), method = method)
   corLinks = melt(corMat)
   
   # construct the links
   tf_source = actLinks[which(actLinks$value > miTh),1]
   tf_target = actLinks[which(actLinks$value > miTh),2]
   tf_links = data.frame(as.matrix(NetMat[NetMat$from %in% tf_source & NetMat$to %in% tf_target, ]), 
                         stringsAsFactors = F)
   tf_links$relation = NULL
   
   # determine the directions
   for (i in 1:nrow(tf_links)){
      source_name = tf_links$from[i]
      target_name = tf_links$to[i]
      tf_links$relation[i] <- ifelse(corMat[source_name, target_name] > 0, 1, 2)
   }
   rownames(tf_links) = NULL
   return(tf_links)
}









#=========================================================#
#==========    racipe model classification   =============# 
#=========================================================# 


########
# _7: Each model is compared with permutations of random data
# Two new functions ClustFunction and ModelPValue
# Assigns clusters based on p values
# should work with all mean, min, z score type assignments as p values are related
# A separate funtion to change the assignment type. Change ClustFunction

HeatmapSimilarity = function(data.reference, data.simulation, cluster.cut = NULL, n.clusters = 3,
                             p.value=0.05, permuted.var, permutations = 1000,
                             cor.method = "spearman", cluster.method = "ward.D2", method = "pvalue", buffer = 0.01, permutation.method = "simulation") {
   #' Calculates the similarity between two gene expression data
   #'
   #'  Comparison is done across columns, i.e., how similar are the columns in the two dataset.
   #' For gene expression data, format data so that gene names are in rows and samples in columns.
   #' @param data.reference The reference data matrix, for example, the experimental gene expression values
   #' @param data.simulation The data matrix to be compared.
   #' @param n.clusters (optional) The number of clusters in which the reference data should be clustered for comparison.
   #' Not needed if cluster.cut is provided
   #' @param p.value (optional) p-value to consider two gene expression sets as belonging to same cluster.
   #' Ward's method with spearman correlation is used to determine if a model belongs to a specific cluster.
   #' @param permuted.var (optional) Similarity scores computed after permutations.
   #' @param cluster.cut (optional) Clsuter numbers assigned to reference data.
   #'  If cluster.cut is missing, hierarchical clustering using /code{ward.D2}
   #'  and /code{distance  = (1-cor(x, method = "spear"))/2} will be used to cluster the reference data.
   #'  @param cluster.method (optional) Clustering method to be used to cluster the experimental data
   #'  @param cor.method (optional) Correlation method. Default method is "spearman". For single cell data, use "kendall"
   #'  @param permutation.method "sample" or "reference"
   #' @return A list containing the KL distance of new cluster distribution from reference data and
   #' the probability of each cluster in the reference and simulated data.
   #'
   #'
   
   
   message("Calculating the similarity index")
   #  n.clusters = 3
   n.models <- dim(data.reference)[2]
   n.models.KO <- dim(data.simulation)[2]
   
   if (missing(permutations)) {
      permutations = 1000
   }
   
   if (missing(cor.method)) {
      cor.method <- "spearman"
   }
   
   ref.cor <- cor((data.reference), method = cor.method)
   
   if (missing(cluster.cut)) {
      if(missing(n.clusters)){
         stop("Please specify the number of clusters using n.clusters or
              cluster assignments using cluster.cut")
      }
      
      # cluster the reference data if the clutering assignments has not been provided.
      distance <- as.dist((1-ref.cor)/2)
      clusters <- hclust(distance, method = cluster.method)
      #plot(clusters)
      cluster.cut <- cutree(clusters, n.clusters)
      
      } else {
         if(!missing(n.clusters)){
            warnings("Neglecting n.clusters. The number of clusters will be determined from cluster.cut.")
         }
         n.clusters <- length(unique(cluster.cut))
      }
   
   # find the variance within each cluster
   #TO DO Will standard deviation be better? shouldn't be with ward method.
   
   ref.cluster.var <- c(rep(0,n.clusters))
   for(j in 1:n.clusters)
   {
      temp.cluster.var <- (((1 - ref.cor[which(cluster.cut==j), which(cluster.cut==j)])/2)^2)
      ref.cluster.var[j] <- ClustFunction(temp.cluster.var[upper.tri(temp.cluster.var, diag = FALSE)])
      temp.cluster.var <- NULL
   }
   
   
   #  cluster.cut <- cluster.cut[1:10]
   #  data.reference <- data.reference[,1:10]
   simulated.ref.cor <- t(cor(data.reference, data.simulation, method = cor.method))
   
   #clusterFreq <- table(CLUSTERCUT)/n_models
   
   if (sum(is.na(simulated.ref.cor)) > 0) {
      message("Error in correlation. Please verify the data")
   }
   
   simulated.cluster.var <- matrix(0, nrow=n.models.KO, ncol = n.clusters)
   
   for(i in 1:n.models.KO){
      for(j in 1:n.clusters)
      {
         temp.cluster.var <- ((1 - simulated.ref.cor[i, which(cluster.cut==j)])/2)^2
         simulated.cluster.var[i,j] <- ClustFunction(temp.cluster.var )
         temp.cluster.var <- NULL
      }
   }
   
   # testing clustering robustness
   # ref.cluster.var <- matrix(0, nrow = n.models, ncol =  n.clusters)
   # for(i in 1:n.models)
   # for(j in 1:n.clusters)
   # {
   #   ref.cluster.var[i,j] <- mean(((1 - ref.cor[i, which(cluster.cut==j)])/2)^2)
   # }
   
   if (method == "variance") {
      simulated.cluster <- matrix(0, nrow =  n.models.KO, ncol = 2)
      simulated.cluster[, 2] <- apply(simulated.cluster.var,1,min)
      # simulated.cluster.allowed <- simulated.cluster.var < ref.cluster.var
      simulated.cluster[, 1] <- apply(simulated.cluster.var,1,which.min)
      simulated.cluster[which(3*ref.cluster.var[simulated.cluster[,1]] < simulated.cluster[, 2]), 1] <- 0
      simulated.cluster <- simulated.cluster[,-2]
      
   }
   #  permutations = 1000
   if(missing(method)) {
      method = "pvalue"
   }
   if (method == "pvalue" ) {
      message("pvalue method")
      if(missing(permuted.var )) {
         if(permutation.method == "reference"){
            permuted.var <- PermutedVar(simulated.ref.cor, cluster.cut, permutations, ref.cluster.var)
            simulated.var.P.value <- SimulatedVarPValue(permuted.var, p.value)
            #rowSums(simulated.cluster.allowed)
            #simulated.cluster.var.sorted <- sort(simulated.cluster.var, index.return = TRUE )
            # simulated.cluster.allowed <- simulated.cluster.var < simulated.var.P.value
            simulated.cluster <- matrix(0, nrow =  n.models.KO, ncol = 2)
            simulated.cluster[, 2] <- apply(simulated.cluster.var,1,min)
            simulated.cluster[, 1] <- apply(simulated.cluster.var,1,which.min)
            simulated.cluster[which(simulated.var.P.value[simulated.cluster[,1]] < simulated.cluster[, 2]), 1] <- 0
            simulated.cluster <- simulated.cluster[,-2]
            
         } else {
            message("simulation permutation")
            
            p.value.mat <- ModelPvalue(data.simulation, data.reference, cluster.cut, permutations,
                                       ref.cluster.var, cor.method, simulated.cluster.var)
            simulated.cluster <- matrix(0, nrow =  n.models.KO, ncol = 2)
            simulated.cluster[, 2] <- apply(p.value.mat,1,min)
            simulated.cluster[, 1] <- apply(p.value.mat,1,which.min)
            simulated.cluster[which(simulated.cluster[,2] > p.value), 1] <- 0
            simulated.cluster <- simulated.cluster[,-2]
            
         }
      }
      
      
   }
   
   similarity <- list()
   cluster.names <- unique(cluster.cut)
   print(c("Original Clusters", cluster.names))
   #cluster.names <- c(cluster.names, "4") #test
   
   simulated.cluster.names <- unique(simulated.cluster)
   print(c("Simulated Clusters", simulated.cluster.names))
   missing.ref.clusters <- setdiff(cluster.names, simulated.cluster.names)
   #print(c("Missing Clusters", missing.ref.clusters))
   bufferEnteriesPerCluster <- max(1,as.integer(buffer*n.models.KO))
   missing.ref.clusters.add <- numeric() #c(rep(0,bufferEnteriesPerCluster*length(missing.ref.clusters)))
   if (length(missing.ref.clusters) > 0) {
      for(i in 1:length(missing.ref.clusters))
      {
         missing.ref.clusters.add <- c(missing.ref.clusters.add, rep(missing.ref.clusters[i],bufferEnteriesPerCluster))
      }
   }
   simulated.cluster.adjusted <- c(simulated.cluster, missing.ref.clusters.add)
   
   
   missing.simulated.clusters <- setdiff(simulated.cluster.names, cluster.names)
   bufferEnteriesPerCluster <- max(1,as.integer(buffer*n.models))
   if(bufferEnteriesPerCluster == 0) bufferEnteriesPerCluster <- 1
   missing.ref.clusters.add <- numeric() #c(rep(0,bufferEnteriesPerCluster*length(missing.ref.clusters)))
   if (length(missing.simulated.clusters) > 0) {
      missing.ref.clusters.add <- c(missing.ref.clusters.add, rep(missing.simulated.clusters,bufferEnteriesPerCluster))
   }
   cluster.cut.adjusted <- c(cluster.cut, missing.ref.clusters.add)
   
   
   ref.cluster.freq <- table(cluster.cut.adjusted)/(length(cluster.cut.adjusted))
   similarity$ref.cluster.freq <- table(cluster.cut)/n.models
   simulated.cluster.freq <- table(simulated.cluster.adjusted)/length(simulated.cluster.adjusted)
   
   similarity$simulated.cluster.freq <- table(simulated.cluster)/n.models.KO
   
   similarity$cluster.similarity <- ref.cluster.freq*log(ref.cluster.freq/simulated.cluster.freq)
   similarity$KL <- sum(similarity$cluster.similarity )
   
   similarity$data.reference <- data.reference
   colnames(similarity$data.reference) <- cluster.cut
   similarity$data.reference <- similarity$data.reference[,order(colnames(similarity$data.reference))]
   
   
   similarity$data.simulation <- data.simulation[,which(simulated.cluster>0)]
   colnames(similarity$data.simulation) <- simulated.cluster[which(simulated.cluster>0)]
   similarity$data.simulation <- similarity$data.simulation[,order(colnames(similarity$data.simulation))]
   ref.sim.cor <- numeric()
   previous.cluster.size <- 0
   ref.sim.cor.ref <- numeric()
   previous.cluster.size.ref <- 0
   
   for(i in 1:(length(unique(colnames(similarity$data.simulation)))))
   {
      temp.ref <- similarity$data.reference[,which(colnames(similarity$data.reference)==i)]
      temp.sim <- similarity$data.simulation[,which(colnames(similarity$data.simulation)==i)]
      
      temp.ref.sim.cor <- cor(temp.ref,temp.sim, method = cor.method)
      ref.sim.cor <- c(ref.sim.cor,previous.cluster.size +
                          sort(colMeans(temp.ref.sim.cor), decreasing = T, index.return = T)$ix)
      previous.cluster.size <- previous.cluster.size + dim(temp.sim)[2]
      
      ref.sim.cor.ref <- c(ref.sim.cor.ref, previous.cluster.size.ref +
                              sort(rowMeans(temp.ref.sim.cor), decreasing = T, index.return = T)$ix)
      previous.cluster.size.ref <- previous.cluster.size.ref + dim(temp.ref)[2]
      
      
   }
   similarity$data.simulation <- similarity$data.simulation[,ref.sim.cor]
   similarity$data.simulation <- cbind(similarity$data.simulation[,ref.sim.cor], data.simulation[,which(simulated.cluster == 0)])
   
   similarity$data.reference <- similarity$data.reference[,ref.sim.cor.ref]
   
   #TO DO : This invovlves repeat calculation of cor--can be optimized
   similarity$simulated.ref.cor <- t(cor(similarity$data.reference, similarity$data.simulation, method = cor.method))
   
   #image(similarity$simulated.ref.cor, col = plot_color)
   return(similarity)
}

#########################################################
# Helper functions
#########################################################

NthMin <- function(x,index) {
   #' Find nth minimum value from a vector
   #' @param x the given unsorted vector
   #' @param index
   #' @return the nth minimum element of the vector
   #'
   return (sort(x, decreasing = FALSE, partial = index)[index])
   
}

#############################################

ClustFunction <- function(x){
   #return (mean(x))
   return (min(x))
}

#############################################

PermutedVar <- function(simulated.ref.cor, cluster.cut, permutations, ref.cluster.var){
   #' Returns the variance array after permutations
   #' @param simulated.ref.cor Correlation matrix of simulated and reference data
   #' @param cluster.cut The original cluster assignments
   #' @param permutations The number of permutations
   #' @return An array of dimension n.models by n.clusters by permutations
   #'
   n.clusters <- length(unique(cluster.cut))
   n.models.KO <- dim(simulated.ref.cor)[1]
   permuted.var <- array(0, c(n.models.KO, n.clusters, permutations))
   for(k in 1:permutations){
      cluster.cut.permuted <- sample(cluster.cut)
      for(j in 1:n.clusters)
      {
         cor.mat <- simulated.ref.cor[,which(cluster.cut.permuted == j)]
         for(i in 1:n.models.KO){
            temp.cluster.var <- (((1-cor.mat[i, ])/2)^2)
            permuted.var[i,j,k] <- (mean(temp.cluster.var)/ (ref.cluster.var[j]))
         }
      }
   }
   return(permuted.var)
}
#############################################

ModelPvalue <- function(data.simulation, data.reference, cluster.cut, permutations,
                        ref.cluster.var, cor.method, simulated.cluster.var){
   #' Returns the variance array after permutations
   #' @param data.simulation
   #' @param data.reference
   #' @param cluster.cut The original cluster assignments
   #' @param permutations The number of permutations
   #' @param ref.cluster.var SD of the clusters
   #' @return An array of dimension n.models by n.clusters by permutations
   #'
   n.clusters <- length(unique(cluster.cut))
   n.models.KO <- dim(data.simulation)[2]
   n.gene <- dim(data.simulation)[1]
   p.value.mat <- matrix(0, nrow = n.models.KO, ncol = n.clusters)
   random.models <-  matrix(rep(seq(1:n.gene),permutations), n.gene, permutations)
   random.models <- apply(random.models,2,sample)
   #random.models <- matrix(0, nrow = n.gene, ncol = permutations)
   #random.models <- t(apply(data.simulation,1,function(x)  sample(x, replace = TRUE, size = permutations)))
   permuted.ref.cor <- matrix(0,nrow = permutations, ncol = dim(data.reference)[2])
   permuted.ref.cor <- cor(random.models, data.reference, method = cor.method)
   for(j in 1:n.clusters)
   {
      dist.mat <- ((1 - permuted.ref.cor[,which(cluster.cut == j)])/2)^2
      temp.vector <- sort(apply(dist.mat,1,ClustFunction))
      for (i in 1:n.models.KO) {
         p.value.mat[i,j] <- (which(abs(temp.vector - simulated.cluster.var[i,j])==min(abs(temp.vector - simulated.cluster.var[i,j])))[1] - 1)/permutations
         #[1] as sometimes which() might satisfy for multiple values
      }
      
   }
   
   return(p.value.mat)
}

#############################################

SimulatedVarPValue <- function(permuted.var, p.value){
   #' Finds the variance corresponding to a given value
   #' @param permuted.var An array containing the distance of clusters for each model for every permutation
   #' @param p.value
   #' @return p-values for each model
   #'
   permutations <- dim(permuted.var)[3]
   n.models.KO <- dim(permuted.var)[1]
   n.clusters <- dim(permuted.var)[2]
   selected.index = as.integer(permutations*p.value)
   
   if(selected.index==0) {stop("Number of permutations is not sufficient to achieve the required p.value.
                               Please increase the permutations")}
   #print(selected.index)
   simulated.var.P.value <- matrix(0, nrow=n.models.KO, ncol = n.clusters)
   
   for (i in 1:n.models.KO) {
      for(j in 1:n.clusters)    {
         simulated.var.P.value[i,j] <- NthMin(permuted.var[i,j,],selected.index)
      }
   }
   
   return(simulated.var.P.value)
   }
#############################################
#' Returns absolute P value
SimulatedPValueAbs <- function(permuted.var, simulated.cluster.var){
   #' Finds the variance corresponding to a given value
   #' @param permuted.var An array containing the distance of clusters for each model for every permutation
   #' @param p.value
   #' @return p-values for each model
   #'
   permutations <- dim(permuted.var)[3]
   n.models.KO <- dim(permuted.var)[1]
   n.clusters <- dim(permuted.var)[2]
   
   simulated.var.P.value <- matrix(0, nrow=n.models.KO, ncol = n.clusters)
   
   for (i in 1:n.models.KO) {
      for(j in 1:n.clusters)    {
         temp.vector <- sort(permuted.var[i,j,])
         
         simulated.var.P.value[i,j] <- which(abs(temp.vector - simulated.cluster.var[i,j])==min(abs(temp.vector - simulated.cluster.var[i,j])))/permutations
      }
   }
   
   return(simulated.var.P.value)
}

