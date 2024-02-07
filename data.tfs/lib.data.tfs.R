TF_Activity.steps = function (tfs, GSDB, eset, DErslt, 
                        with_weight = TRUE, 
                        if_module = FALSE, 
                        ind = 1/5){ 
   if (is(eset, "ExpressionSet")){
      data = exprs(eset)
   }else{data = eset}
   table = DErslt$table
   DEweights = data.frame(p = table$padj)
   rownames(DEweights) = rownames(table)
   tfs = tfs[tfs %in% rownames(data)]
   all_activity = all_list = NULL
   
   results.inter <- list() 
   
   print('in TF_Activity function')
   #source('./lib.data.tfs.R')
   library(Hmisc)
   for (tf in tfs) {
      comms = intersect(GSDB[[tf]], rownames(table)) 
      results.inter[['comm']] <- comms 
      
      g_inds = which(rownames(data) %in% comms)
      tmp_data = data[g_inds, ]
      tmp_data = rem_data(tmp_data)
      results.inter[['tmp_data']] <- tmp_data
      
      # TF Activity 1
      #==============
      
      # Construct correlation matrix using target expressions
      #------------------------------------------------------
      cor_mat = rcorr(as.matrix(t(tmp_data)), type = "spearman")[["r"]]   
      dimnames(cor_mat) = list(rownames(tmp_data), rownames(tmp_data))   
      results.inter[['cor_mat']] <- cor_mat


      # Construct column sum vector from correlation matrix  
      #---------------------------------------------------- 
      tmp_names = list(rownames(tmp_data), "sign")
      tmp_dist = cor_mat
      col_sums = apply(tmp_dist, 1, sum) 
      
      results.inter[['col_sums']] <- col_sums
      
      # Construct mod matrix by scaling correlation matrix 
      #---------------------------------------------------
      mod_mat = tmp_dist - sapply(col_sums, function(x) x *
                                     col_sums/sum(tmp_dist))  
      results.inter[['mod_mat']] <- mod_mat 

      # Calculate eigen values and eigen vectors of mod matrix 
      #-------------------------------------------------------
      ev = eigen(mod_mat) 
      tmp_sign = data.frame(sign(ev$vectors[, 1]))
      dimnames(tmp_sign) = tmp_names
      ng_mat = cor_mat[tmp_sign$sign == -1, tmp_sign$sign ==-1]
      ps_mat = cor_mat[tmp_sign$sign == 1, tmp_sign$sign == 1]
      test_mats = list(ng_mat = ng_mat, ps_mat = ps_mat) 
      results.inter[['ev']] <- ev
      results.inter[['tmp_sign']] <- tmp_sign
      results.inter[['test_mats']] <- test_mats
   
      # Retain the targets that are within the mean distance  
      # from the center for each of +ve and -ve matrix
      #----------------------------------------------------
      gs_list1 = list()
      test_center = numeric()
      for (i in 1:2) {
         test_mat = test_mats[[i]]
         nrow_mat = nrow(test_mat)
         if (is.null(nrow_mat)) {
            next
         }
         else if (nrow_mat <= 2) {
            gs_list1[[i]] = rownames(test_mat)
         }
         else {
            tmp_kmean = kmeans(test_mat, 1)
            centers = tmp_kmean$centers
            test_center = c(test_center, centers)
            test_dim = length(centers)
            distances = apply(test_mat, 1, function(x) cor(as.numeric(x),
                                                           as.numeric(centers), method = "spearman"))
            m = mean(distances)
            names_left = names(distances)[distances >= m]
            gs_list1[[i]] = names_left
         }
      }
      results.inter[['gs_list1']] <- gs_list1

      
      # Calculate TF activity based on these remaining targets
      #----------------------------------------------------------------------
      tf_exprs = as.numeric(data[tf, ])
      gs_remain = unlist(gs_list1)
      tmp_rslt = cal_activity(gs_remain, tmp_data, tmp_sign,
                              ind, with_weight, DEweights, tf_exprs)
      tmp_activity = tmp_rslt$activity
      tmp_sign = tmp_rslt$sign 
      
      results.inter[['TF_actvity_1']] <- list('act'=tmp_activity, 'sign'= tmp_sign)
      
      
      # TF Activity 2
      #==============
      # TF activity 2: calculate TF activity based on target correlations
      #------------------------------------------------------------------ 
      gs_list2 = list()
      gs_remain2 = c()
      cors = apply(tmp_data, 1, function(x) cor(x, tf_exprs,
                                                method = "spearman"))
      pos_cors = cors[cors > 0]
      neg_cors = cors[cors < 0]
      tmp_sign2 = data.frame(sign(cors))
      m_pos = mean(pos_cors)
      m_neg = mean(neg_cors)
      gs_list2[[1]] = names(pos_cors[pos_cors >= m_pos])
      gs_list2[[2]] = names(neg_cors[neg_cors <= m_neg]) 
      
      results.inter[['gs_list_2']] <- gs_list2
      
      gs_remain = unlist(gs_list2)
      tmp_rslt = cal_activity(gs_remain, tmp_data, tmp_sign2,
                              ind, with_weight, DEweights, tf_exprs)
      tmp_activity2 = tmp_rslt$activity
      tmp_sign2 = tmp_rslt$sign
      
      # Calculate spearman correlations between target expressions 
      # and TF activity 1
      #-----------------------------------------------------------
      cor1 = apply(tmp_data, 1, function(x) {
         tmpc = cor(x, tmp_activity, method = "spearman")
         return(tmpc)
      })
      
      # Calculate spearman correlations between target expressions 
      # and TF activity 2
      #-----------------------------------------------------------
      cor2 = apply(tmp_data, 1, function(x) {
         tmpc = cor(x, tmp_activity2, method = "spearman")
         return(tmpc)
      }) 
      
      results.inter[['cor1']] <- cor1
      results.inter[['cor2']] <- cor2
 
      
      # Decide which activity to report based on comparative 
      # values of the two mean correlations   
      #----------------------------------------------------- 
      if (if_module == T){
         all_list[[tf]] = tmp_sign
         all_activity = rbind(all_activity, tmp_activity)
      }else{
         mean1 = mean(abs(cor1))
         mean2 = mean(abs(cor2))
         if (mean1 < mean2) {
            tmp_sign = tmp_sign2
            tmp_activity = tmp_activity2
         }
         all_list[[tf]] = tmp_sign
         all_activity = rbind(all_activity, tmp_activity)
      } 
      results.inter[['tmp_sign']] <- tmp_sign 
      results.inter[['tmp_activity']] <- tmp_activity
      return(results.inter)
   }
   # dimnames(all_activity) = list(tfs, colnames(data))
   # return(list(all_list = all_list, all_activities = all_activity[complete.cases(all_activity),
   #                                                               ]))
}

select.balanced.tfs.by_pval <- function (tfs.df, 
                                         GSDB, 
                                         eset, 
                                         de.results, 
                                         COLNAME.PREFIX.CONTROL, 
                                         COLNAME.PREFIX.TREATMENT, 
                                         PVALUE.CUTOFF){
   
   # Select balanced TFs for CTRL-IDH1 group
   #----------------------------------------- 
   library(NetAct) 
   
   # Calculate activities for CTRL-IDH1 
   #------------------------------------------------ 
   a = TF_Activity(tfs = as.character(tfs.df$tf),
                   GSDB=GSDB,  
                   eset=eset,
                   DErslt=de.results
   )
   
   # Compare Average expression for CONTROL and TREATMENT 
   #-----------------------------------------------------
   new_activity <- a$all_activities
   class(new_activity)
   dim(new_activity)
   tfs.df 
   # TFS lost in the activity calculations: 
   #length(as.character(tfs.ctrl.idh1.df$tf)) - dim(new_activity)[1]
   length(as.character(tfs.df$tf)) - dim(new_activity)[1]
   
   # Select TOP TFS upregulated in CONTROL and TREATMENT 
   #----------------------------------------------------
   # Retain rows that correspond to TFs returned by TF_Activity 
   gs = rownames(new_activity)
   data <- data[gs, ]
   dim(data)
   
   # Normalize data 
   data.norm <- row_norm(data)
   
   # Separate CONTROL and TREATMENT data 
   data.norm.ctrl <- as.data.frame(data.norm) %>% dplyr::select(starts_with(COLNAME.PREFIX.CONTROL))
   data.norm.trtmt <- as.data.frame(data.norm) %>% dplyr::select(starts_with(COLNAME.PREFIX.TREATMENT))
   
   dim(data.norm.ctrl)
   dim(data.norm.trtmt)
   
   # Calculate average expressions of genes across all CONTROLs:  
   avg_expr.ctrl <- rowMeans(data.norm.ctrl)
   # Calculate average expressions of genes across all TREATMENTs:
   avg_expr.trtmt <- rowMeans(data.norm.trtmt)
   
   length(avg_expr.ctrl)
   length(avg_expr.trtmt)
   
   head(avg_expr.ctrl)
   head(avg_expr.trtmt)
   
   # Find genes upregulated in CONTROL:
   UP.genes.CTRL <- row.names(data.norm)[avg_expr.ctrl > avg_expr.trtmt]
   length(UP.genes.CTRL)
   
   # Find genes upregulated in TREATMENT: 
   UP.genes.TRTMT <- row.names(data.norm)[avg_expr.ctrl < avg_expr.trtmt]
   length(UP.genes.TRTMT) 
   
   head(UP.genes.CTRL)
   head(UP.genes.TRTMT)
   
   # Retain TFs that satisfy p-value criteria
   #------------------------------------------
   tfs.pval.df <- tfs.df[which(tfs.df$pvals<=PVALUE.CUTOFF), ]
   
   # Collect UP TFs from UP genes for CTRL and TRTMT
   #------------------------------------------------
   UP.TFs.CTRL.pval <- tfs.pval.df[which(tfs.pval.df$tf %in% UP.genes.CTRL), ]
   UP.TFs.TRTMT.pval <- tfs.pval.df[which(tfs.pval.df$tf %in% UP.genes.TRTMT), ]
   
   # Retain TFs that satisfy p-value criteria
   #------------------------------------------
   TF.count.CTRL <- nrow(UP.TFs.CTRL.pval)
   TF.count.TRTMT <- nrow(UP.TFs.TRTMT.pval)
   
   cat(TF.count.CTRL)
   cat(TF.count.TRTMT)
   
   if(TF.count.CTRL<TF.count.TRTMT){
      UP.TFs.CTRL <- UP.TFs.CTRL.pval[1:TF.count.CTRL, ]
      UP.TFs.TRTMT <- UP.TFs.TRTMT.pval[1:TF.count.CTRL, ]
   }else if(TF.count.CTRL>TF.count.TRTMT){
      UP.TFs.CTRL <- UP.TFs.CTRL.pval[1:TF.count.TRTMT, ]
      UP.TFs.TRTMT <- UP.TFs.TRTMT.pval[1:TF.count.TRTMT, ]
   }else{
      UP.TFs.CTRL <- UP.TFs.CTRL.pval
      UP.TFs.TRTMT <- UP.TFs.TRTMT.pval
   } 
   
   nrow(UP.TFs.CTRL)
   nrow(UP.TFs.TRTMT)
   
   balanced.tfs <- list()
   balanced.tfs[['CTRL']] <- UP.TFs.CTRL
   balanced.tfs[['TRTMT']] <- UP.TFs.TRTMT   
   
   return(balanced.tfs)
}


select.balanced.tfs <- function(tfs.df, 
                                GSDB, 
                                eset, 
                                de.results, 
                                COLNAME.PREFIX.CONTROL, 
                                COLNAME.PREFIX.TREATMENT, 
                                NO_TFS.CONTROL, 
                                NO_TFS.TREATMENT){
   library(NetAct) 
   
   # Calculate activities for CTRL-IDH1 
   #------------------------------------------------ 
   a = TF_Activity(tfs = as.character(tfs.df$tf),
                   GSDB=GSDB, #GSDB=hDB,
                   eset=eset,
                   DErslt=de.results
   )
   
   
   # Compare Average expression for CONTROL and TREATMENT 
   #-----------------------------------------------------
   new_activity <- a$all_activities
   class(new_activity)
   dim(new_activity)
   tfs.df 
   # TFS lost in the activity calculations: 
   #length(as.character(tfs.ctrl.idh1.df$tf)) - dim(new_activity)[1]
   length(as.character(tfs.df$tf)) - dim(new_activity)[1]
   
   # Select TOP TFS upregulated in CONTROL and TREATMENT 
   #----------------------------------------------------
   # Retain rows that correspond to TFs returned by TF_Activity 
   gs = rownames(new_activity)
   data <- data[gs, ]
   dim(data)
   
   # Normalize data 
   data.norm <- row_norm(data)
   
   # Separate CONTROL and TREATMENT data 
   data.norm.ctrl <- as.data.frame(data.norm) %>% select(starts_with(COLNAME.PREFIX.CONTROL))
   data.norm.trtmt <- as.data.frame(data.norm) %>% select(starts_with(COLNAME.PREFIX.TREATMENT))
   
   dim(data.norm.ctrl)
   dim(data.norm.trtmt)
   
   # Calculate average expressions for CONTROL and TREATMENT
   avg_expr.ctrl <- rowMeans(data.norm.ctrl)
   avg_expr.trtmt <- rowMeans(data.norm.trtmt)
   
   length(avg_expr.ctrl)
   length(avg_expr.trtmt)
   
   head(avg_expr.ctrl)
   head(avg_expr.trtmt)
   
   # Find TFs upregulated in CONTROL
   tfs.UP.ctrl <- row.names(data.norm)[avg_expr.ctrl > avg_expr.trtmt]
   length(tfs.UP.ctrl)
   
   # Find TFs upregulated in TREATMENT 
   tfs.UP.trtmt <- row.names(data.norm)[avg_expr.ctrl < avg_expr.trtmt]
   length(tfs.UP.trtmt) 
   
   head(tfs.UP.ctrl)
   head(tfs.UP.trtmt)
   
   # Select equal number of TFs from TFs upregulated in CONTROL and TREATMENT
   top_tfs.UP.ctrl <- head(tfs.UP.ctrl, NO_TFS.CONTROL)
   top_tfs.UP.trtmt <- head(tfs.UP.trtmt, NO_TFS.TREATMENT)
   length(top_tfs.UP.ctrl)
   length(top_tfs.UP.trtmt)
   
   # top_tfs <- union(top_tfs.UP.ctrl, top_tfs.UP.trtmt)
   # length(top_tfs)
   # length(unique(top_tfs))   
   
   balanced.tfs <- list()
   balanced.tfs[['CTRL']] <- top_tfs.UP.ctrl 
   balanced.tfs[['TRTMT']] <- top_tfs.UP.trtmt 
   return(balanced.tfs)
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
   #celltypes = c(rep("CTRL", 11), rep("IDH1", 10), rep("IDH2", 10))
   celltypes = c(rep("CTRL", 11), rep("IDH", 9))
   
   phenoData = new("AnnotatedDataFrame", data = data.frame(celltype = celltypes,  
                                                           samples = colnames(edata))) 
   rownames(phenoData) = colnames(edata)
   
   # Construct an instance of ExpressionSet 
   #---------------------------------------
   eset = ExpressionSet(assayData = as.matrix(edata), phenoData = phenoData) 
   
   return(eset)
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
