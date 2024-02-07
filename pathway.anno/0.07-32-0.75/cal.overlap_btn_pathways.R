remove(list = ls()) 

NO_PATHWAYS <- 10 
SIG_LEVEL <- 0.05

source('../functions.R')


outdir <- './pathway.overlaps/'
dir.create(outdir)

# Load circuit DEGS
#------------------
degs.circuit <- read.table(file = './degs/DEGs_in_circuit.txt', header = T)
degs.circuit <- as.character(degs.circuit$x)
length(degs.circuit)

nTOTAL <- length(degs.circuit)

# Load circuit DEGS
#------------------
degs.pathways <- read.table(file = './enrichr.degs_in_circuit/KEGG_2019_Human_table.txt', 
                            header = T, sep = '\t') 

degs.pathways <- degs.pathways[1:NO_PATHWAYS, ]
rownames(degs.pathways) <- degs.pathways$Term


# Find overlap between DEGs in cellcycle and cellular senescence
#----------------------------------------------------------------
geneset.cellcycle <- as.character(degs.pathways[2,"Genes"]) 
geneset.cellularSen <- as.character(degs.pathways[8,"Genes"])

tmp1 <- strsplit(geneset.cellcycle, split=';', fixed = TRUE)[[1]] 

tmp2 <- strsplit(geneset.cellularSen, split=';', fixed = TRUE)[[1]]
length(tmp2)

overlapGenes <- intersect(tmp1, tmp2)
length(overlapGenes)

length(geneset.cellcycle)

cat(overlapGenes)


# Drop unwanted pathways: Pathways in cancer, Epstein-Barr virus infection, and Hepatitis B
#--------------------------------------------------------------------------------------
row.names(degs.pathways)
#pathways_to_remove <- c("Pathways in cancer", "Epstein-Barr virus infection", "Hepatitis B", "Cellular senescence")
pathways_to_remove <- c("Pathways in cancer", "Epstein-Barr virus infection", "Hepatitis B")
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

 
# Fisher's Exact Test (FET): find pathways that are very similar
#--------------------------------------------------------------
pathway.names <- rownames(degs.pathways)
n = length(pathway.names)
stat_pair.wise_pathway.enrichment <- as.data.frame(matrix(nrow = n*(n-1)/2, ncol = 7))
colnames(stat_pair.wise_pathway.enrichment) <- c('pathway_1', 'pathway_2', 'pvalue', 
                                                 'no_tfs.pathway_1', 'no_tfs.pathway_2', 
                                                 'no_overlapped_tfs', 'overlapped_tfs')

count <- 1
for(pno_1 in 1:(length(pathway.names)-1)){ 
   print(pno_1) 
   pname_1 <- rownames(degs.pathways)[pno_1]
   print(pname_1)
   degs.pathway_1 <- as.character(degs.pathways[pname_1, "Genes"]) 
   degs.pathway_1 <- unlist(strsplit(degs.pathway_1, "\\;")) 
   for(pno_2 in (pno_1+1):length(pathway.names)){ 
      print(pno_2) 
      pname_2 <- rownames(degs.pathways)[pno_2]
      print(pname_2)
      degs.pathway_2 <- as.character(degs.pathways[pname_2, "Genes"]) 
      degs.pathway_2 <- unlist(strsplit(degs.pathway_2, "\\;"))   
      
      # perform Fisher's exact test:
      ftest <- FET_1not2(glist1=degs.pathway_1, glist2=degs.pathway_2, ntot=nTOTAL) 
      #ftest <- FET_1not2(glist1=degs.pathway_2, glist2=degs.pathway_1, ntot=nTOTAL) 

      tfs.overlapped <- intersect(degs.pathway_1, degs.pathway_2)
      
      stat_pair.wise_pathway.enrichment[count, ] <- c(pname_1, pname_2, format(ftest$pval, digits = 2), 
                                                      length(degs.pathway_1), length(degs.pathway_2), 
                                                      length(tfs.overlapped), 
                                                      paste(tfs.overlapped, sep = ' ', collapse = ';'))
      count <- count + 1
      #break()
   } 
   #break()
}


fname.out <- paste(outdir, 'stat_pair.wise_pathway.enrichment.csv', sep = '')
write.csv(stat_pair.wise_pathway.enrichment, file = fname.out, row.names = F, quote = F)


# Select the pathways 
stat_pair.wise_pathway.enrichment$pvalue <- as.numeric(stat_pair.wise_pathway.enrichment$pvalue)

# Set H0 and Ha
# H0: the two pathways are enriched 
#     i.e. the association is NOT random 
# Ha: two pathways are different

# test for H0
# if pvalue <= SIG_LEVEL, we fail to reject H0 ==> not enough evidence to support that 
# the two pathways are different
stat_H0 <- stat_pair.wise_pathway.enrichment[stat_pair.wise_pathway.enrichment$pvalue<=SIG_LEVEL,]  

# test Ha
# if pvalue > SIG_LEVEL, we reject H0 in favor of Ha ==> enough evidence to support that 
# the two pathways are different
stat_Ha <- stat_pair.wise_pathway.enrichment[stat_pair.wise_pathway.enrichment$pvalue>SIG_LEVEL,]


fname.out <- paste(outdir, 'stat_pair.wise_H0.csv', sep = '')
write.csv(stat_H0, file = fname.out, row.names = F, quote = F)

fname.out <- paste(outdir, 'stat_pair.wise_Ha.csv', sep = '')
write.csv(stat_Ha, file = fname.out, row.names = F, quote = F)


# Obtain the pathways from Ha
#-----------------------------
pathways.Ha <- unique(unique(stat_Ha$pathway_1), unique(stat_Ha$pathway_2))
cat(pathways.Ha)
length(pathways.Ha) 
NULL

