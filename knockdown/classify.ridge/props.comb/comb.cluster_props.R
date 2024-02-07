remove(list = ls()) 

datadir <- './data/'
dir.create(datadir)

# Load cluser proportions for single and double knockdown 
#------------------------------------------------------- 
fname.in <- '../props.skd/data/cluster_props.csv'
cluster_props.skd <- read.csv(file = fname.in, row.names = 1) 

fname.in <- '../props.dkd/data/cluster_props.csv'
cluster_props.dkd  <- read.csv(file = fname.in, row.names = 1) 

# Combine both single and double knockdown cluster props  
#-------------------------------------------------------  
cluster_props.dkd.new <-  cluster_props.dkd[!rownames(cluster_props.dkd) %in% c("Untreated"),] 
cluster_props.comb <- rbind(cluster_props.skd, cluster_props.dkd.new)
 
fname.out <- paste(datadir, 'cluster_props.csv', sep = '')
write.csv(cluster_props.comb, file = fname.out, quote = F)
