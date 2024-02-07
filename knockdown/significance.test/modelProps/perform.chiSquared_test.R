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


# Perform Chi-squared tests
#--------------------------
stats.TF <- as.data.frame(matrix(nrow = (nrow(cluster_props.sorted)-1), 
                                 ncol = (ncol(cluster_props.sorted)+3)) )
colnames(stats.TF) <- c('TF', colnames(cluster_props.sorted), 'pvalue', 'nlogpvalue')
rownames(stats.TF) <- setdiff(rownames(cluster_props.sorted), c("Untreated")) 

for(tf in setdiff(rownames(cluster_props.sorted), c("Untreated"))){
   print(tf)
   chisq <- chisq.test(as.numeric(cluster_props.sorted[tf,])*100,
                       p=as.numeric(cluster_props.sorted["Untreated",]))
   stats.TF[tf, ] <- c(tf, as.numeric(cluster_props.sorted[tf,]), 
                       chisq$p.value, -log10(chisq$p.value))  
   #break()
}

# Convert the values to numeric 
#-------------------------------
stats.TF$AML <- as.numeric(stats.TF$AML)
stats.TF$`Hybrid+Untreated` <- as.numeric(stats.TF$`Hybrid+Untreated`)
stats.TF$pvalue <- as.numeric(stats.TF$pvalue)
stats.TF$nlogpvalue <- as.numeric(stats.TF$nlogpvalue)

# Save data
#----------
fname.out <- paste(outdir , 'stats.Chi-squared_test.csv', sep = '')
write.csv(stats.TF, file = fname.out, quote = F) 


