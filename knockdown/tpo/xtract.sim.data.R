#!/usr/bin/env Rscript

remove(list = ls()) 

circuit_tpo <- read.table('./circuit.tpo', header=T)
class(circuit_tpo)
head(circuit_tpo)

length(unique(sort(c(circuit_tpo$Source, circuit_tpo$Target)) ))

# extract TFs in the network 
#---------------------------
TFs_in_circuit <- sort(union(circuit_tpo$Source, circuit_tpo$Target))
length(TFs_in_circuit)

fname.out <- 'TFs-in-circuit.txt'
write.table(TFs_in_circuit, file = fname.out, sep = '\t', quote = F, row.names = F)
