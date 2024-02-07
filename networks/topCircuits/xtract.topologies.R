#!/usr/bin/env Rscript

remove(list = ls())

#circuit_idx <- '0.07-32-0.75'
circuit_names <- c('circuit_1', 'circuit_2', 'circuit_3', 'circuit_4') 

circuit_indices <- c('0.07-32-0.75', '0.07-32-0', '0.09-32-0.85', '0.16-28-0.85') 
# names(circuit_indices) <- circuit_names
names(circuit_names) <- circuit_indices

library(sRACIPE) 

# outdir <- './tpo/'
# dir.create(outdir)
outdir <- './'

# Extact circuits
#----------------- 
nodes.all <- c()
tpo.list <- list()
nodes.byCircuit <- list()
for(circuit_idx in names(circuit_names)){
  racipe <- readRDS(file = paste('../circuits.sim/circuit_simulated_', 
                                 circuit_idx, '.rds', sep = '') )
  circuit_tpo <- sracipeCircuit(racipe)
  fname.out <- paste(outdir, 'circuit-', circuit_idx, '.tpo' ,sep = '')
  write.table(circuit_tpo, file = fname.out, sep = '\t', quote = F, row.names = F) 
  
  tpo.list[[circuit_names[circuit_idx]]] <- circuit_tpo
  
  TFs_in_circuit <- union(circuit_tpo$Source, circuit_tpo$Target)
  nodes.byCircuit[[circuit_names[circuit_idx]]] <- TFs_in_circuit 
  nodes.all <- c(nodes.all, TFs_in_circuit)
}

names(tpo.list)
names(nodes.byCircuit) 
nodes.all <- unique(nodes.all) # 59


# find common nodes in each circuit pair
#---------------------------------------
common.nodes.df <- as.data.frame(matrix(nrow = 6, ncol = 4))
colnames(common.nodes.df) <- c('circuit.pair', 'num_common.nodes', 'nodeCount.pair', 'circuit.idx.pair')
count.pair <- 1
for(n1 in 1:(length(circuit_names)-1)){
  circuit.1 <- as.character(circuit_names[n1])
  for(n2 in (n1+1):(length(circuit_names))){
    circuit.2 <- as.character(circuit_names[n2]) 
    common.tfs <- intersect(nodes.byCircuit[[circuit.1]], nodes.byCircuit[[circuit.2]]) 
    common.nodes.df[count.pair, ] <-  c(paste(circuit.1, '.', circuit.2, sep = ''), 
                                        length(common.tfs), 
                                        paste(length(nodes.byCircuit[[circuit.1]]), '_', length(nodes.byCircuit[[circuit.2]]), sep = ''),
                                        paste(names(circuit_names[n1]), '_', names(circuit_names[n2]), sep = ''))
    count.pair <- count.pair + 1
    #break()
  }
  #break()
}

write.csv(common.nodes.df, file = paste(outdir, 'common.nodes.csv', sep = ''), quote = F, row.names = F)


# get number of nodes in each circuit
#------------------------------------
nodes.df <- as.data.frame(matrix(nrow = 4, ncol = 3))
colnames(nodes.df) <- c('circuit.name', 'num_nodes', 'circuit.idx')
count.circuit <- 1
for(n1 in 1:(length(circuit_names))){
  circuit.name <- as.character(circuit_names[n1]) 
  print(circuit.name) 
  nodes.df[count.circuit, ] <- c(circuit.name,  length(nodes.byCircuit[[circuit.name]]),  names(circuit_names[n1]))
  count.circuit <- count.circuit + 1

}
write.csv(nodes.df, file = paste(outdir, 'nodes.perCircuit.csv', sep = ''), quote = F, row.names = F)


# Assign status to each node of the largest circuit (circuit 2)
#-------------------------------------------------------------- 
# summary.status.df <- as.data.frame(matrix(nrow = 4, ncol = 3))
# colnames(summary.status.df) <- c('circuit.name', 'circuit.index', 'num_nodes', '')

nodes.all <- unique(nodes.all)
nodes.all <- sort(nodes.all)
node.status.df <- as.data.frame(matrix(nrow = length(nodes.all), ncol = 2))
colnames(node.status.df) <- c('tf', 'status')
node.status.df$tf <- nodes.all
node.status.df$status <- 4

# assign status 1 to nodes of circuit 4
node.status.df[node.status.df$tf %in% nodes.byCircuit$circuit_4, 2] <- 1
sum(node.status.df$status==1) # 22

# ckt 4 and ckt 3
tmp1 <- node.status.df
new.nodes <- setdiff(nodes.byCircuit$circuit_3, nodes.byCircuit$circuit_4) 
length(new.nodes) # 8
node.status.df[node.status.df$tf %in% new.nodes, 2] <- 2
sum(node.status.df$status==2)


# ckt 3 and ckt 1
tmp2 <- node.status.df
new.nodes <- setdiff(nodes.byCircuit$circuit_1, nodes.byCircuit$circuit_3) 
length(new.nodes)
new.nodes <- setdiff(new.nodes, nodes.byCircuit$circuit_4) 

node.status.df[node.status.df$tf %in% new.nodes, 2] <- 3
sum(node.status.df$status==3) # 23


# ckt 1 and ckt 2
tmp3 <- node.status.df
new.nodes <- setdiff(nodes.byCircuit$circuit_2, nodes.byCircuit$circuit_1) 
length(new.nodes) # 6
new.nodes <- setdiff(new.nodes, nodes.byCircuit$circuit_1) 
length(new.nodes) # 6

#node.status.df[node.status.df$tf %in% new.nodes, 2] <- 4
sum(node.status.df$status==4) # 6

sum(node.status.df$status==1) # 22
sum(node.status.df$status==2) # 8
sum(node.status.df$status==3) # 23
sum(node.status.df$status==4) # 6

sum(nodes.byCircuit$circuit_2 %in% 'TP63') # 0

# remove one node that is not in circuit 2
dim(node.status.df) # 59 x 2 
sum(node.status.df$tf %in% nodes.byCircuit$circuit_2) 
node.status.df <- node.status.df[node.status.df$tf %in% nodes.byCircuit$circuit_2, ] 
dim(node.status.df) # 58 x 2 - TP63 is gone

write.csv(node.status.df, file = paste0(outdir, 'node.status.csv'), quote = F, row.names = F)



