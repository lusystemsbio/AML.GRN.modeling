
degs.cur <- read.csv(file = './degs/DEGs.txt', header = F) 
degs.phase.37 <- read.csv(file = '/Users/kateba/research/aml.idh/phase.37/networks/pathway.anno/degs/DEGs.txt', 
                          header = F)

degs.cur <- as.character(unlist(degs.cur))
degs.phase.37 <- as.character(unlist(degs.phase.37))


common.degs <- intersect(degs.cur, degs.phase.37)
length(degs.cur) # 8497
length(degs.phase.37) # 8507
length(common.degs) # 8413 