






# Constants
#-------------
COLNAME.PREFIX.CONTROL <- 'CTRL' 
COLNAME.PREFIX.TREATMENT <- 'IDH'
#PVALUE.CUTOFF <- 0.10

# Load expression file
#------------------------- 
# library(Biobase)
# fname.eset <- '../data.tfs/eset.brain_array.rda'
# load(fname.eset) 
# class(eset.brain_array) 
# data <- exprs(eset.brain_array) 
# dim(data)

#tfs.ordered <- as.character(coreTFs$Netact)
tfs.ordered <- as.character(coreTFs$COMB)
hgs <- targetDB.list$`0.05`


upGeneCount <- cal.UPgeneCount(tfs=tfs.ordered, 
                               GSDB=hgs,
                               eset=data,
                               de.results=de.results,
                               COLNAME.PREFIX.CONTROL=COLNAME.PREFIX.CONTROL,
                               COLNAME.PREFIX.TREATMENT=COLNAME.PREFIX.TREATMENT
                               ) 

# set parameters
#----------------
tfs=tfs.ordered 
GSDB=hgs 
eset=data 
# de.results=de.results
# COLNAME.PREFIX.CONTROL=COLNAME.PREFIX.CONTROL
# COLNAME.PREFIX.TREATMENT=COLNAME.PREFIX.TREATMENT
# PVALUE.CUTOFF=PVALUE.CUTOFF


cal.UPgeneCount <- function(tfs, 
                            GSDB, 
                            eset, 
                            de.results, 
                            COLNAME.PREFIX.CONTROL, 
                            COLNAME.PREFIX.TREATMENT){ 
   # Select balanced TFs for CTRL-IDH1 group
   #----------------------------------------- 
   library(NetAct) 
   print(tfs)
   
   # Calculate activities for CTRL-IDH1 
   #------------------------------------------------ 
   a = TF_Activity(tfs=tfs,  
                   GSDB=GSDB,  
                   eset=eset,
                   DErslt=de.results
   ) 
   
   # Compare Average expression for CONTROL and TREATMENT 
   #-----------------------------------------------------
   new_activity <- a$all_activities
   class(new_activity)
   dim(new_activity)
   tfs 
   
   # TFS lost in the activity calculations: 
   length(as.character(tfs)) - dim(new_activity)[1]
   
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
   
   UPgeneCount <- list() 
   UPgeneCount[['CTRL']] <- length(UP.genes.CTRL)
   UPgeneCount[['TRTMT']] <- length(UP.genes.TRTMT)
   return(UPgeneCount)
}   


#==============================================
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
   a = TF_Activity(tfs = tfs.df, #as.character(tfs.df$tf),
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
   #length(as.character(tfs.df$tf)) - dim(new_activity)[1]
   length(as.character(tfs.df)) - dim(new_activity)[1]
   
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


