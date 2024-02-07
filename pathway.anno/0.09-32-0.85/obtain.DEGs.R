#!/usr/bin/env Rscript

remove(list = ls()) 

outdir <- './degs/'
dir.create(outdir)

# DE results
#------------
fname.de.results <- '../../data.tfs/de.results.rda'
load(fname.de.results) # loads object de.results 
class(de.results)

# extract DE genes: 
DEGs <- de.results$degs
class(DEGs) 
rm(de.results)
length(DEGs)
write.table(DEGs, file = paste(outdir, 'DEGs.txt', sep = ''), 
            col.names = F, row.names = F, quote = F)

