#!/usr/bin/env Rscript

remove(list = ls()) 

library(sRACIPE)
library(NetAct)
# source('../heatmapSimilarity.updated.R')
# source('../functions.R')
library(gplots)

# Constants 
#----------
circuit_idx <- '0.09-32-0.85'
NO_MODELS <- 10000

# input directory for hS objects obtained from cal.metrics.sim_circuits.R
hS.dir <- '../../networks/circuits.hS/' 

outdir <- './data/'
dir.create(outdir)
figdir <- './figs/'
dir.create(figdir)

# Explore one selected circuit
#-----------------------------
racipe <- readRDS(file = paste('../../networks/circuits.sim/circuit_simulated_', 
                               circuit_idx, '.rds', sep = '') )
names(racipe)
length(names(racipe))

# Plot circuit
#-------------
#sRACIPE::sracipePlotCircuit(racipe, plotToFile = F) 
sRACIPE:: sracipePlotCircuit(racipe, plotToFile = T) 

class(racipe)
names(racipe)
assays(racipe)

circuit_tpo <- sracipeCircuit(racipe)
length(unique(sort(c(circuit_tpo$Source, circuit_tpo$Target)) ))

# Save the circuit
#----------------- 
#outdir <- './data/'
#fname.out <- paste(outdir, 'circuit-', circuit_idx, '.tpo' ,sep = '')
#write.table(circuit_tpo, file = fname.out, sep = '\t', quote = F, row.names = F)

# extract TFs in the network 
#---------------------------
TFs_in_circuit <- union(circuit_tpo$Source, circuit_tpo$Target)
length(TFs_in_circuit)

fname.out <- paste(outdir, 'TFs-', circuit_idx, '.txt' ,sep = '')
write.table(TFs_in_circuit, file = fname.out, sep = '\t', quote = F, row.names = F)

# Plot simulation data
#---------------------
sRACIPE::sracipePlotData(racipe, plotToFile = T)
# sRACIPE::sracipePlotData(racipe, plotToFile = F) 

# Load heatmap similarity object
#-------------------------------
fname.hS <- paste(hS.dir, 'hS_', circuit_idx, '.rds', sep = '')
fname.hS
hS <- readRDS(file = fname.hS)
hS$AvgDist
hS$simulated.cluster.freq
hS$simulated.cluster.freq[2] + hS$simulated.cluster.freq[3]
hS$KL

# Help for heatmap.2 graphical params: https://www.biostars.org/p/312405/ 
# par(oma=c(1,1,1,1)); # bottom, left, top, right
# heatmap.2(hS$dataReference, trace = 'none') 

# Heatmap for reference data
#---------------------------
hS.ref <- hS$dataReference
WIDTH <- 8
HEIGHT <- 8
#figname <- paste(figdir, 'heatmap.refData-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
#pdf(file = figname, paper = 'special', width = WIDTH, height = HEIGHT) 

figname <- paste(figdir, 'heatmap.refData-.pdf', sep = '') 
pdf(file = figname, paper = 'special') 
par(mfrow = c(1, 1)) 
# par(mar=c(12.1,2.1,2.0,2.1)) # bottom, left
par(oma=c(2,1,1,1)); # bottom, left, top, right
heatmap.2(hS$dataReference, trace = 'none') 
dev.off() 

# Plot simulated.refCor 
#---------------------
# hS$simulated.refCor: pearson correlation between data.REF and data.sim 
simulated.refCor <- hS$simulated.refCor
dim(simulated.refCor)

WIDTH <- 8
HEIGHT <- 6
figname <- paste(figdir, 'heatmap.simulated.refCor-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')
pdf(file = figname, paper = 'special', width = WIDTH, height = HEIGHT)
graphics::image(hS$simulated.refCor) 
#graphics::image(t(hS$simulated.refCor))
dev.off()

# Heatmap for simulation data
#----------------------------
col2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
          "#D55E00", "#CC79A7") 
plotColor <- c("#5E4FA2", "#4F61AA", "#4173B3", "#3386BC", "#4198B6",
               "#51ABAE", "#62BEA6", "#77C8A4", "#8ED1A4", "#A4DAA4",
               "#B8E2A1", "#CBEA9D", "#DEF199", "#EAF69F", "#F2FAAC",
               "#FAFDB8", "#FEFAB6", "#FEF0A5", "#FEE695", "#FDD985",
               "#FDC978", "#FDB96A", "#FCA75E", "#F99254", "#F67D4A",
               "#F26943", "#E85A47", "#DE4B4B", "#D33C4E", "#C1284A",
               "#AF1446", "#9E0142")
simClusters <- as.character(col2[(1+hS$simClusters)]) 

WIDTH <- 8
HEIGHT <- 6
# figname <- paste(figdir, 'heatmap.simulated-', WIDTH, 'x', HEIGHT,'.pdf', sep = '')
# pdf(file = figname, paper = 'special', width = WIDTH, height = HEIGHT) 

figname <- paste(figdir, 'heatmap.simulated.pdf', sep = '')
pdf(file = figname, paper = 'special')
gplots::heatmap.2(hS$dataSimulation, trace = "none",
                  dendrogram = "none", Colv=FALSE, col = plotColor,
                  ColSideColors = simClusters,
                  main = "Simulated Data",
                  distfun=function(x) as.dist(1-cor(t(x), method = "s")))
dev.off()
