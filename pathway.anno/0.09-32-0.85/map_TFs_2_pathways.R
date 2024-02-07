#!/usr/bin/env Rscript

remove(list = ls()) 

NO_PATHWAYS <- 12 # 10
SIGN.LEVEL <- 0.1 #0.05 #0.01

# outdirectory
outdir <- './tfs.annotated/'
dir.create(outdir)

source('./functions.R') 

# Load circuit DEGS
#------------------
degs.circuit <- read.table(file = './degs/DEGs_in_circuit.txt', header = F)
degs.circuit <- as.character(degs.circuit$V1)
length(degs.circuit)

nTOTAL <- length(degs.circuit)
nTOTAL

# Load pathways and involved DEGS
#---------------------------------
degs.pathways <- read.table(file = './enrichr.degs_in_circuit/KEGG_2019_Human_table.txt', 
                            header = T, sep = '\t') 


degs.pathways <- degs.pathways[1:NO_PATHWAYS, ]
rownames(degs.pathways) <- degs.pathways$Term

# drop two pathways: Pathways in cancer, Epstein-Barr virus infection, and Hepatitis B
#--------------------------------------------------------------------------------------
#pathways_to_remove <- c("Pathways in cancer", "Epstein-Barr virus infection", "Hepatitis B") 
pathways_to_remove <- c("Pathways in cancer", "Epstein-Barr virus infection", "Hepatitis B", 
                        "Measles", "Human papillomavirus infection") 

degs.pathways <- degs.pathways[!rownames(degs.pathways) %in% pathways_to_remove, ]

NO_PATHWAYS <- NO_PATHWAYS - length(pathways_to_remove)
NO_PATHWAYS 

degs.pathways$Genes[1]
length(degs.pathways$Genes[1])

degs.pathway <- degs.pathways$Genes[1]

class(degs.pathway)
degs.pathway <- as.character(degs.pathway)
degs.pathway  <- unlist(strsplit(degs.pathway, "\\;"))

length(degs.pathway)
head(degs.pathway)


# Load TOP circuit 
#-----------------
#circuit_metrics.sim <- read.csv(file = '../networks/results/summary.circuits.sim.sorted.csv', row.names = 1)
#circuit_idx <- rownames(circuit_metrics.sim)[1]
circuit_idx <- "0.09-32-0.85"
circuit_idx

fname.hS <- paste('../../networks/circuits.hS/hS_', circuit_idx, '.rds', sep = '')
hS <- readRDS(file = fname.hS)
class(hS)
names(hS)
data.sim <- hS$dataSimulation
class(data.sim)
data.ref <- hS$dataReference
dim(data.ref)

TFs_in_circuit <- rownames(data.ref)
length(TFs_in_circuit) 

# Load TF DEGs
#-------------
# all target DBs at each feature ratio cutoff 
targetDB.list <- readRDS(file = '../../databases/targetDB.list.rds')
length(names(targetDB.list)) 

# select the target DB specific to the feature ratio used 
# to create the current circuit: 
regDB.feature.ratio <- feature_ratio_cutoff <- strsplit(circuit_idx, '-', 2)[[1]][1]  
regDB.feature.ratio
targetDB <- targetDB.list[[regDB.feature.ratio]]
length(names(targetDB))

# Limit target DB (both TFs and their targets) to the DEGs related to circuit
#-----------------------------------------------------------------------
# limit TFs of the target DB  to the DEGs related to the circuit:
targetDB <- targetDB[TFs_in_circuit] 
length(names(targetDB))

# limit targets of the target DB  to the DEGs related to the circuit:
for(tf in names(targetDB)){
   targetDB[[tf]] <- intersect(degs.circuit, targetDB[[tf]])   
   #break()
}


# Fisher's Exact Test (FET) for all selected pathways 
#----------------------------------------------------
tf.sig.table <- as.data.frame(matrix(nrow = dim(degs.pathways)[1], ncol = length(targetDB)))
rownames(tf.sig.table) <- degs.pathways$Term
colnames(tf.sig.table) <- names(targetDB) 

for(pname in rownames(degs.pathways)){
   #print(pname) 
   degs.pathway <- as.character(degs.pathways[pname, "Genes"]) 
   #degs.pathway  <- unlist(strsplit(degs.pathway, "\\;")) 
   pathway_DEGS <- unlist(strsplit(degs.pathway, "\\;")) 
   for(tf in names(targetDB)){
      TF_DEGS <- c(tf, targetDB[[tf]])
      ftest <- FET_1not2(glist1=TF_DEGS, glist2=pathway_DEGS , ntot=nTOTAL) 
      #ftest <- FET_1not2(glist1=pathway_DEGS, glist2=TF_DEGS, ntot=nTOTAL) 
      tf.sig.table[pname, tf] <- ftest$pval
   }
   #break()
}

dim(tf.sig.table)

# Select significant genes for each pathway
#------------------------------------------
tf.sig.table.t <- t(tf.sig.table)
tf.sig.table.t <- tf.sig.table.t[order(rownames(tf.sig.table.t)), ]

fname.out <- paste(outdir, 'tf.sig.table.csv', sep = '')
write.csv(format(tf.sig.table.t, digits = 2), file = fname.out, quote = F, row.names = T)

# Select the top pathways for each TF
#------------------------------------
sig_patways_byTF <- as.data.frame(matrix(nrow = dim(tf.sig.table.t)[1], ncol = (NO_PATHWAYS*2+1)))
colnames(sig_patways_byTF) <- c('tf', paste('pathway_', rep(1:NO_PATHWAYS, times=1), sep = ''),
                                        paste('pathway_', rep(1:NO_PATHWAYS, times=1), '.pvalue', sep = '') )

rownames(sig_patways_byTF) <- rownames(tf.sig.table.t)

for(tf in rownames(tf.sig.table.t)){
   #print(tf)
   pathways.sig <- tf.sig.table.t[tf,]
   pathways.sig.sorted <- sort(pathways.sig)
   sig_patways_byTF[tf, ] <- c(tf, names(pathways.sig.sorted), as.numeric(pathways.sig.sorted))
   #break
}


TFsWithAnnoted_pathways <- as.data.frame(matrix(nrow = length(rownames(sig_patways_byTF)), ncol = 3))
colnames(TFsWithAnnoted_pathways) <- c('tf','count','pathways') 
rownames(TFsWithAnnoted_pathways) <- rownames(sig_patways_byTF)
count.outer <- 1
count.inner <- 1
for(tf in rownames(sig_patways_byTF)){
   #print(tf)
   patways.sele <- c()
   for(pathway.no in 1:NO_PATHWAYS){
      colname.pathway <- paste('pathway_', pathway.no, sep = '')
      colname.pvalue <- paste('pathway_', pathway.no, '.pvalue', sep = '')
      #if(as.numeric(sig_patways_byTF[tf, colname.pvalue])<=SIGN.LEVEL){
         pathway.name <- sig_patways_byTF[tf, colname.pathway] 
         pathway.pvalue <- sig_patways_byTF[tf, colname.pvalue]
         pathway.pvalue <- format(as.numeric(pathway.pvalue), digits = 1)
         # pathway.pvalue.pair <- paste(sig_patways_byTF[tf, colname.pathway], '(', sig_patways_byTF[tf, colname.pvalue], ')' ,sep = '') 
         pathway.pvalue.pair <- paste(pathway.name, '(', pathway.pvalue, ')' ,sep = '') 
         patways.sele <- c(patways.sele, pathway.pvalue.pair)
      #}
      count.inner <- count.inner + 1
   } 
   patways.sele.collapsed <- paste(patways.sele, collapse = ';', sep = '')
   TFsWithAnnoted_pathways[tf, ] <- c(tf, length(patways.sele), patways.sele.collapsed)
   count.outer <- count.outer + 1
   #break()
}

TFsWithAnnoted_pathways <- TFsWithAnnoted_pathways[order(TFsWithAnnoted_pathways$tf), ]

fname.out <- paste(outdir, 'TFsWithAnnoted_pathways.csv', sep = '') 
fname.out
write.csv(TFsWithAnnoted_pathways, file = fname.out, quote = F, row.names = F)



#============= part 2 ====================================# 
# Create TF to pathway id mapping 
#-------------------------------- 
pathway_names <- sort(c(as.character(degs.pathways$Term))) 
# pathway_names <- c(pathway_names, 'NA') 
# length(pathway_names)
# pathway_ids <- c(1, 2, 2, 3, 4, 5, 6, 7)
# names(pathway_ids) <- pathway_names

# remove MAPK - because this pathway is not significant (pvalue 0.1)
# for any TF
pathway_names <- setdiff(pathway_names, pathway_names[5]) 
pathway_names <- c(pathway_names, 'NA') 
length(pathway_names)
pathway_ids <- c(1, 2, 2, 3, 4, 5, 6)
names(pathway_ids) <- pathway_names


# Save the TF to pathway id mapping 
#----------------------------------- 
pathway_vs_id <- cbind(names(pathway_ids), as.character(pathway_ids))
colnames(pathway_vs_id) <- c('pathway', 'pathway_id')

fname.out <- paste(outdir, 'pathway_vs_id', '.csv', sep = '') 
fname.out
write.csv(pathway_vs_id, file = fname.out, quote = F, row.names = F)


# Annotate each TF in the circuit with only one pathway or None  
# Annotate each TF with most most significant pathway
#--------------------------------------------------------------
TFsWith_singleAnnotedPathway <- as.data.frame(matrix(nrow = length(TFs_in_circuit), ncol = 3))
colnames(TFsWith_singleAnnotedPathway) <- c('tf', 'pathway', 'pathway_id')
rownames(TFsWith_singleAnnotedPathway) <- TFs_in_circuit

for(tf in TFs_in_circuit){
   #print(tf) 
   pathways.cur <- unlist(strsplit(TFsWithAnnoted_pathways[tf, "pathways"], "\\;"))
   pathway.top <- unlist(strsplit(pathways.cur[1], "\\("))[1] 
   pvalue.tmp <-  unlist(strsplit(pathways.cur[1], "\\("))[2] 
   pvalue.top <- as.numeric(unlist(strsplit(pvalue.tmp, "\\)")))
   if(!is.na(print(pathway.top)) && pvalue.top <= SIGN.LEVEL){
      TFsWith_singleAnnotedPathway[tf, ] <- c(tf, names(pathway_ids[pathway.top]), 
                                                  as.character(pathway_ids[pathway.top]))
   }else {
      TFsWith_singleAnnotedPathway[tf, ] <- c(tf, names(pathway_ids['NA']), 
                                              as.character(pathway_ids['NA']))
   }
   #break()
}

TFsWith_singleAnnotedPathway <- TFsWith_singleAnnotedPathway[order(TFsWith_singleAnnotedPathway$tf), ]

# number of times each pathway is associated with any TF
sum(TFsWith_singleAnnotedPathway$pathway_id==1) # 2 AMPK 
sum(TFsWith_singleAnnotedPathway$pathway_id==2) # 12 Cell cycle/Cellular senescence  
sum(TFsWith_singleAnnotedPathway$pathway_id==3) # 4 JAK-STAT
sum(TFsWith_singleAnnotedPathway$pathway_id==4) # 7 p53 
sum(TFsWith_singleAnnotedPathway$pathway_id==5) # 10 AMPK
sum(TFsWith_singleAnnotedPathway$pathway_id==6) # 3 NA


# Save the TF to pathway mapping
#-------------------------------
fname.out <- paste(outdir, 'TFsWith_singleAnnotedPathway.csv', sep = '') 
fname.out
write.csv(TFsWith_singleAnnotedPathway, file = fname.out, quote = F, row.names = F)

