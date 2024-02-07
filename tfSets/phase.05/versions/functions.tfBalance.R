selectTOPtfs.byMethod <- function(TFs.combined, targetDB, top.TFs.count){ 
   # TFs.combined=TFs.combined
   # targetDB=targetDB.cur
   # top.TFs.count=top.TFs.count
   
   # select TFs from Netact
   tfs.df=TFs.combined[['Netact']][,c(1, 6)] # select columns: tfs and qvalue  
   tfs.df <-  tfs.df[order(tfs.df[,2], decreasing = FALSE), ] # order them by qvalue - smallest to largest
   tfs.Netact <- selectTOPtfsSigVal(tfs.df=tfs.df, targetDB=targetDB, MIN_TARGETS=8)  
   
   # select TFs from MARINa
   tfs.df=TFs.combined[['MARINa']][,c(1, 4)] # select columns: tfs and FDR
   tfs.df <-  tfs.df[order(tfs.df[,2], decreasing = FALSE), ] # order them by FDR - smallest to largest
   tfs.MARINa <- selectTOPtfsSigVal(tfs.df=tfs.df, targetDB=targetDB, MIN_TARGETS=8)  
   
   # select TFs from RI
   tfs.df=TFs.combined[['RI']][,c(1, 2)] # select columns: tfs and aggregate error
   tfs.df <-  tfs.df[order(tfs.df[,2], decreasing = TRUE), ] # order them by FDR - largest to smallest
   tfs.RI <- selectTOPtfsSigVal(tfs.df=tfs.df, targetDB=targetDB, MIN_TARGETS=8)  
   
   tfs.byMethod <- list()
   tfs.byMethod[["Netact"]] <- tfs.Netact 
   tfs.byMethod[["MARINa"]] <- tfs.MARINa
   tfs.byMethod[["RI"]] <- tfs.RI
   tfs.byMethod[['COMB']] <- unique(sort(union(union(as.character(tfs.Netact), 
                                               as.character(tfs.MARINa)), 
                                               as.character(tfs.RI))))
   return(tfs.byMethod)
}


#===========
selectTOPtfsSigVal <- function(tfs.df, targetDB, MIN_TARGETS=8){ 
   # MIN_TARGETS <- 8
   # method.name <- names(TFs.combined)[1] # Netact
   # tfs.df <- TFs.combined[[method.name]] 
   # tfs.df <- tfs.df[,c(1, 6)] # tested for Netact
   
   rownames(tfs.df) <- tfs.df$tf
   # Retain the tfs with targets >= MIN_TARGETS - filtering out TFs with too few targets
   tfs.allowed <- as.data.frame(matrix(nrow = 0, ncol = ncol(tfs.df))) 
   colnames(tfs.allowed) <- colnames(tfs.df) 
   for(tf in rownames(tfs.df)){
      if(length(targetDB[[tf]]) >= MIN_TARGETS) { 
         tfs.allowed <- rbind(tfs.allowed, tfs.df[tf, ]) 
      } 
   } 
   colnames(tfs.allowed) <- colnames(tfs.df) 
   
   # select top TFs 
   if(dim(tfs.allowed)[1]>=top.TFs.count){
      tfs.top <- as.character(tfs.allowed[1:top.TFs.count, 1])
      names(tfs.top) <- tfs.allowed[1:top.TFs.count, 2] # name the tfs with sig. values
   }else{
      tfs.top <- as.character(tfs.allowed[, 1])
      names(tfs.top) <- tfs.allowed[,2]
   }  
   return(tfs.top)
}


#==========================
# combine TOP TFs from each method from the pool of TFs  
# from each method that satisfy minimum number of targets criteria
selectTOPtfs <- function(TFs.combined, targetDB, top.TFs.count){
   MIN_TARGETS <- 8
   tfs.sele <- c()
   for(method.name in names(TFs.combined)){
      tfs.cur <- as.character(TFs.combined[[method.name]][,1]) 
      tfs.ordered <- c()
      for(x in tfs.cur) {
         if(length(targetDB[[x]]) >= MIN_TARGETS) {
            tfs.ordered <- c(tfs.ordered, x)
         }
      }
      tfs.sele <- c(tfs.sele, 
                    tfs.ordered[which(complete.cases(tfs.ordered[1:top.TFs.count]))])
   }
   return(tfs.sele)
}



plotStat.comb <- function(stats.df, colors, figname, ylimit){
   colors <- c('black', 'magenta','blue', 'red') 
   names(colors) <- c("TFs", "ACTUAL", "CTRL", "TRTMT") 
   WIDTH <- 10
   HEIGHT <- 8 #16 #8 
   pdf(file = figname, width=WIDTH, height=HEIGHT, paper='special')
   plot(stats.df$no_tfs*3, col=colors["TFs"], pch = 19, type = 'n',
        xlim = c(1, 90), ylim = ylimit,
        xlab = 'Index', ylab = 'Number of TFs')
   #points(stats.df$no_tfs_actual, col=colors["ACTUAL"])
   points(stats.df$CTRL, col=colors["CTRL"])
   points(stats.df$TRTMT, col=colors["TRTMT"])
   dev.off()   
}


plotStat <- function(stats.df, colors, figname, ylimit){
   colors <- c('black', 'magenta','blue', 'red') 
   names(colors) <- c("TFs", "ACTUAL", "CTRL", "TRTMT")  
   WIDTH <- 10
   HEIGHT <- 8 
   pdf(file = figname, width=WIDTH, height=HEIGHT, paper='special')
   plot(stats.df$no_tfs, col=colors["TFs"], pch = 19, type = 'n',
        xlim = c(1, 90), ylim = ylimit, 
        xlab = 'Index', ylab = 'Number of TFs')
   #points(stats.df$no_tfs_actual, col=colors["ACTUAL"])
   points(stats.df$CTRL, col=colors["CTRL"])
   points(stats.df$TRTMT, col=colors["TRTMT"]) 
   dev.off()   
}


cal.UPgeneCount <- function(tfs, 
                            GSDB, 
                            eset, 
                            de.results, 
                            COLNAME.PREFIX.CONTROL, 
                            COLNAME.PREFIX.TREATMENT){ 
   # Select balanced TFs for CTRL-IDH1 group
   #----------------------------------------- 
   library(NetAct) 
   #print(tfs)
   
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





# Select the TOP n number of TFs from the combined TFs 
# from each of three methods: NetAct, MARINa, and RI 
# Save TFs from individual method and combined
sel.topTFs_forChecking <- function(TFs.combined, 
                                   NO_TOP_TFS_FROM_EACH_METHOD = 30){ 
   
   # TFs.combined=TFs.combined
   # NO_TOP_TFS_FROM_EACH_METHOD=top.TFs.count
   
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
   
   tfs.sele <- list('Netact' = top.tfs.netact, 
                    'MARINa' = top.tfs.MARINa, 
                    'RI' = top.tfs.RI,
                    'COMB' = tfs.COM)
   
   return(tfs.sele)
}
