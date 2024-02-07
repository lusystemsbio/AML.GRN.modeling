#!/usr/bin/env Rscript
remove(list = ls()) 

circuit_metrics.sim <- read.csv(file = './results/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)
dim(circuit_metrics.sim)

figdir <- './figs.distr/'
dir.create(figdir)
fname_fig <- paste0(figdir, 'networksize-vs-accuracy.pdf')
pdf(file=fname_fig, width=6, height=8, paper = "special", onefile = TRUE) 
par(mfrow=c(3,1))
par(mar=c(5.0, 5.5, 1.5, 2)) # bottom, left, top, right
plot(circuit_metrics.sim$Nodes, circuit_metrics.sim$Accuracy, 
     xlab='Nodes', ylab='')

plot(circuit_metrics.sim$Interactions, circuit_metrics.sim$Accuracy, 
     xlab='Interactions', ylab='')

plot(circuit_metrics.sim$PosInt, circuit_metrics.sim$Accuracy, 
     xlab='Positive Interactions', ylab='')
dev.off()


class(circuit_metrics.sim$TopTFs)




# Bar plot - accuracy 
#---------
data.bp <- cbind(rownames(circuit_metrics.sim),
                 circuit_metrics.sim$Nodes,
                 circuit_metrics.sim$Accuracy) 

colnames(data.bp) <- c('circuit_idx',  'nodes', 'accuracy') 
data.bp <- as.data.frame(data.bp)

data.bp$nodes <- sprintf("%03d", as.integer(data.bp$nodes))
data.bp$nodes <- as.factor(data.bp$nodes)

data.bp$accuracy <- as.numeric(data.bp$accuracy)


library(reshape2)
data_long <- melt(data.bp, id.vars = c('circuit_idx', 'nodes' )) 
colnames(data_long)

HEIGHT <- 8

p.nodes <- ggplot(data_long, aes(x=nodes, y=value)) +
           geom_boxplot() + 
           ylim(c(0,1)) + 
           theme(axis.title.y=element_blank(), 
                 axis.text.x = element_text(angle = 90, vjust=0.5))
p.nodes

figname <- paste(figdir, 'nodes_vs_accuracy.pdf', sep = '')
ggsave(filename =  figname, width = 14, height = HEIGHT) 




# Bar plot - flexibility
#---------
data.bp <- cbind(rownames(circuit_metrics.sim),
                 circuit_metrics.sim$Nodes,
                 circuit_metrics.sim$flexibility) 

colnames(data.bp) <- c('circuit_idx',  'nodes', 'flexibility') 
data.bp <- as.data.frame(data.bp)

data.bp$nodes <- sprintf("%03d", as.integer(data.bp$nodes))
data.bp$nodes <- as.factor(data.bp$nodes)

data.bp$flexibility <- as.numeric(data.bp$flexibility)


library(reshape2)
data_long <- melt(data.bp, id.vars = c('circuit_idx', 'nodes' )) 
colnames(data_long)

HEIGHT <- 8

p.nodes <- ggplot(data_long, aes(x=nodes, y=value)) +
  geom_boxplot() + 
  #ylim(c(0,1)) + 
  theme(axis.title.y=element_blank(), 
        axis.text.x = element_text(angle = 90, vjust=0.5))
p.nodes

figname <- paste(figdir, 'nodes_vs_flexibility.pdf', sep = '')
ggsave(filename =  figname, width = 14, height = HEIGHT) 


