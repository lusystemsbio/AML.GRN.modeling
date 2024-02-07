tmp <- cor(hS$dataReference,
           hS$dataSimulation,
           method = 'spearman')
dim(hS$dataReference)
dim(hS$dataSimulation)
dim(tmp)
image(tmp)
image(t(tmp))

tmp <- t(cor(hS$dataReference,
             hS$dataSimulation,
             method = 'spearman'))
image(tmp)



