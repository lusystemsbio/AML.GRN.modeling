#!/usr/bin/env Rscript
# Chi-squared test for goodness of fit

remove(list = ls()) 

datadir <- './data/'
outdir <- './data/'

fname.in <- paste(datadir , 'cluster_props.sorted.csv', sep = '')
mydata <- read.csv(file = fname.in, row.names = 1)

# AML vs (Hybrid, Untreated)
#---------------------------
cluster_props.sorted <- cbind(mydata[, 2], mydata[,3] + mydata[,4])  
colnames(cluster_props.sorted) <- c('AML', 'Hybrid+Untreated')
rownames(cluster_props.sorted) <- rownames(mydata)


# Normalize the values 
#---------------------
for(tf in rownames(cluster_props.sorted)){
   cluster_props.sorted[tf, ] <- as.numeric(cluster_props.sorted[tf,])/sum(as.numeric(cluster_props.sorted[tf,]))
}

# Calculate fold change
#----------------------
# fold change of AML w.r.t to the hybrid and Untreated combined
foldChange <- cluster_props.sorted[, 1]/cluster_props.sorted[, 2]

foldChange.df <- cbind(rownames(cluster_props.sorted), foldChange)
colnames(foldChange.df) <- c('tf', 'foldChange')

# Save data
#----------
fname.out <- paste(outdir , 'foldChange.csv', sep = '')
write.csv(format(foldChange.df, digits = 4), file = fname.out, quote = F, row.names = F)  
