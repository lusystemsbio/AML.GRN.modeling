#!/usr/bin/env Rscript

# Plot statistics for only 10K models
remove(list = ls()) 

NO_TOP_CIRCUITS <- 5
NO_SAMPLES_PER_CIRCUIT <- 10

outdir <- paste('./figs.10K.TOP', NO_TOP_CIRCUITS, '/', sep = '')
dir.create(outdir)

library(reshape2) 
library(gridExtra)
library(ggplot2) 

# Obtain metrics of top 5 circuits across samples for each sample size level
#---------------------------------------------------------------------------
top.circuit_id.list <- list()  
metric.top_circuits.list <- list()

source('./functions.R')

fname.input <- paste("./samples/comb.metrics.SSIZE.10K", '/rankedCircuits.acc.flex/', 
                     'metric.comb.sorted.csv', sep = '') 
print(fname.input)
metric.comb <- read.csv(file = fname.input, row.names = 1)  
metric.top_circuits <- metric.comb[1:NO_TOP_CIRCUITS, 1:NO_SAMPLES_PER_CIRCUIT]  

top.circuit_id.list[['10K']] <- rownames(metric.top_circuits)
metric.top_circuits.list[['10K']] <- metric.top_circuits 

# Obtain the mapping between circuit id and circuit name
#--------------------------------------------------------
mapping_circuit_ids <- create_circuit_id_2_name_mapping(top.circuit_id.list)

# save the mapping
fname.out <- paste(outdir, 'mappings.circuit_id_vs_circuit_name.csv' ,sep = '')
write.csv(mapping_circuit_ids, file = fname.out, row.names = F)

# Obtain blox plots for metrics of top circuits across sample sizes
#-----------------------------------------------------------------
DOT.SIZE <- 0.7 
SIZE.TICK_LABEL <- 12

dotplots_by_samplSize <- list() 

sampleSize <- names(metric.top_circuits.list)[1]

print(sampleSize) 
metric.top_circuits <- metric.top_circuits.list[[sampleSize]]
circuit_ids <- rownames(metric.top_circuits)
circuit_ids
circuit_names <- sapply(circuit_ids, function(circuit_id) 
   as.character(mapping_circuit_ids$circuit_name[mapping_circuit_ids$circuit_id %in% circuit_id])
)
rownames(metric.top_circuits) <- circuit_names


metric.top_circuits.tmp <-as.data.frame(cbind(rownames(metric.top_circuits), 
                                              metric.top_circuits)) 
colnames(metric.top_circuits.tmp) <- c('circuits', 
                                       colnames(metric.top_circuits))
class(metric.top_circuits.tmp$circuits) 
metric.top_circuits.tmp$circuits <- factor(metric.top_circuits.tmp$circuits, 
                                           levels = metric.top_circuits.tmp$circuits) 

metric.top_circuits.m <- melt(metric.top_circuits.tmp, id=c("circuits"))

class(metric.top_circuits.m$value)
#metric.top_circuits.m$value <- as.numeric(metric.top_circuits.m$value)

p1 <- ggplot(metric.top_circuits.m, aes(x = circuits, y = value)) + 
   #geom_boxplot() 
   geom_dotplot(binaxis = 'y', stackdir = 'center', 
                dotsize=DOT.SIZE) + 
   ylim(5, 14) #ylim(0, 25)

p2 <- p1 + theme(text = element_text(size=SIZE.TICK_LABEL), 
                 axis.text.x = element_text(angle = 90, vjust = 0.5), 
                 axis.title.x = element_blank(), axis.title.y = element_blank()) 

# add mean and sd with a cross bar
# p3 <- p2 + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
#                         geom="crossbar", width=0.5) 
# add mean and sd with a pointrange
p3 <-  p2 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                         geom="pointrange", color="red", size=0.4)   

p3 <- p3 + ggtitle(sampleSize) + theme(plot.title = element_text(hjust = 0.5)) 

# turn of y axis tick labels for 2nd and subsequent panels
#if (count.panel> 1) p3 <- p3 + theme(axis.text.y = element_blank())

dotplots_by_samplSize[[sampleSize]] <- p3

p3 

# Save the circuits 
#-------------------
# g <- grid.arrange(grobs=dotplots_by_samplSize, nrow=2)

WIDTH <- 2
HEIGHT <- 3
fname.out <- paste(outdir, 'dotplot_scores-TOP', 
                   NO_TOP_CIRCUITS, 'circuits-', WIDTH, 'x', HEIGHT, '.pdf', 
                   sep = '')
ggsave(filename = fname.out, dotplots_by_samplSize$`10K`, 
       width = WIDTH, height = HEIGHT)

WIDTH <- 6
HEIGHT <- 6
fname.out <- paste(outdir, 'dotplot_scores-TOP', 
                   NO_TOP_CIRCUITS, 'circuits-', WIDTH, 'x', HEIGHT, '.pdf', 
                   sep = '')
ggsave(filename = fname.out, dotplots_by_samplSize$`10K`, 
       width = WIDTH, height = HEIGHT) 


#============================================
# t-test to determine whether 1st circuit is significantly different from 
# the subsequent circuits
#-------------------------------------------------------------------------#
SIG.LEVEL_t.test <- 0.05 
stats.t.test <- NULL
#cnames.stats.t.test <- c('circuits', 'pvalue', 'log10Pvalue') 
cnames.stats.t.test <- c('circuits', 'log10Pvalue') 
dotplots_t.test_by_samplSize <- list()
count.panel <- 1

for (sampleSize in names(metric.top_circuits.list)){
   print(sampleSize) 
   metric.top_circuits <- metric.top_circuits.list[[sampleSize]]
   circuit_ids <- rownames(metric.top_circuits)
   circuit_ids
   circuit_names <- sapply(circuit_ids, function(circuit_id) 
      as.character(mapping_circuit_ids$circuit_name[mapping_circuit_ids$circuit_id %in% circuit_id])
   )
   rownames(metric.top_circuits) <- circuit_names 
   
   scores.circuit_1 <- as.numeric(metric.top_circuits[1,])
   
   # perform t test
   #----------------
   stats.t.test <- NULL
   for(circuit_name in rownames(metric.top_circuits)[2:NO_TOP_CIRCUITS]){
      print(circuit_name) 
      scores.circuit_x <- as.numeric(metric.top_circuits[circuit_name,]) 
      
      t.test.res <- t.test(scores.circuit_1, scores.circuit_x, var.equal = TRUE)
      #stats.t.test <- rbind(stats.t.test, c(sampleSize, circuit_name, t.test.res$p.value)) 
      stats.t.test <- rbind(stats.t.test, c(circuit_name, 
                                            #t.test.res$p.value, 
                                            -log10(t.test.res$p.value)))  
      colnames(stats.t.test) <- cnames.stats.t.test 
   }
   
   # Create plots for pvalues from t test
   stats.t.test <- as.data.frame(stats.t.test)
   #stats.t.test$pvalue <- as.numeric(as.character(stats.t.test$pvalue))
   stats.t.test$log10Pvalue <- as.numeric(as.character(stats.t.test$log10Pvalue))
   stats.t.test$circuits <- factor(stats.t.test$circuits, 
                                   levels = stats.t.test$circuits) 
   
   p1 <- ggplot(stats.t.test, aes(x = circuits, y = log10Pvalue)) + 
      geom_dotplot(binaxis = 'y', stackdir = 'center', 
                   dotsize=DOT.SIZE) + 
      ylim(0, 3.0) #ylim(0, 8)
   p2 <- p1 + theme(text = element_text(size=SIZE.TICK_LABEL), 
                    axis.text.x = element_text(angle = 90, vjust = 0.5), 
                    axis.title.x = element_blank(), axis.title.y = element_blank())
   p3 <- p2 + geom_hline(yintercept = -log10(SIG.LEVEL_t.test), linetype="longdash")
   p4 <- p3 + 
      ggtitle(sampleSize) + theme(plot.title = element_text(hjust = 0.5)) + 
      scale_x_discrete(breaks=stats.t.test$circuits,
                       labels= paste('circuit_1 vs ', 
                                     as.character(stats.t.test$circuits), sep = ''))
   
   # turn of y axis tick labels for 2nd and subsequent panels
   if (count.panel> 1) p4 <- p4 + theme(axis.text.y = element_blank())
   
   dotplots_t.test_by_samplSize[[sampleSize]] <-  p4
   count.panel <- count.panel + 1
   #break()
}

p4 

# arrange the list of figures on a grid
g.test <- grid.arrange(grobs=dotplots_t.test_by_samplSize, nrow=1)

WIDTH <- 2
HEIGHT <- 3

fname.out <- paste(outdir, 'dotplot_log10Pval-TOP', NO_TOP_CIRCUITS, 'circuits-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
ggsave(filename = fname.out, g.test, width = WIDTH, height = HEIGHT)

# Combine both plots: scores and p-values
#------------------------------------------
dotplots.comb <- c(dotplots_by_samplSize, dotplots_t.test_by_samplSize)
g.both <- grid.arrange(grobs=dotplots.comb, nrow=2)
WIDTH <- 6
HEIGHT <- 8
fname.out <- paste(outdir, 'dotplots_scores_and_pvalues-TOP', NO_TOP_CIRCUITS, 'circuits-', WIDTH, 'x', HEIGHT, '.pdf', sep = '')
ggsave(filename = fname.out, g.both, width = WIDTH, height = HEIGHT)
 

WIDTH <- 3
HEIGHT <- 8
fname.out <- paste(outdir, 'dotplots_scores_and_pvalues-TOP', NO_TOP_CIRCUITS, 'circuits-', WIDTH, 'x', HEIGHT, '.pdf', sep = '')
ggsave(filename = fname.out, g.both, width = WIDTH, height = HEIGHT)
