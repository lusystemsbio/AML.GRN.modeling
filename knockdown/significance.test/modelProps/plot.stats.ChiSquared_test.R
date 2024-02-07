
remove(list = ls()) 

datadir <- './data/'

figdir <- './figs.sig.test/'
dir.create(figdir)

# Load data
#----------
fname.data <- paste(datadir , 'stats.Chi-squared_test.csv', sep = '')
stats.TF <- read.csv(file = fname.data, row.names = 1) 

# Prepare Chi-squared test data for plotting
#-------------------------------------------
stats.TF$TF <- factor(stats.TF$TF, levels = stats.TF$TF)

# Plot the log p values 
#-----------------------
SIG.LEVEL <- 0.05
YLIMIT <- c(0.0, 2.0)
OFFSET.VERTICAL <- 0.08

print(-log10(SIG.LEVEL))

#plot(stats.TF$nlogpvalue, ylim=YLIMIT, xlab='', ylab='', xaxt='n')
#abline(h=-log10(SIG.LEVEL))  
#axis(side = 1, at=1:dim(stats.TF)[1], las=2, 
#     labels = rownames(stats.TF), tck = -0.03)
#title(ylab = '-log10 p-value')


library(gridExtra)
library(ggplot2) 
DOT.SIZE <- 0.7 
SIZE.TICK_LABEL <- 12
YLIMIT <- c(0.0,  2.71) #2.65) # max(stats.TF$nlogpvalue)  

plot_title <- paste('-log (pvalue)', 
                    sep = '  ') 

p1 <- ggplot(stats.TF, aes(x = TF, y = nlogpvalue)) + 
   geom_dotplot(binaxis = 'y', stackdir = 'center', 
                dotsize=DOT.SIZE) + 
   ylim(YLIMIT) + 
   coord_flip() 
p1

# p2 <- p1 + theme(text = element_text(size=SIZE.TICK_LABEL), 
#                  #axis.text.x = element_text(angle = 90, vjust = 0.5), 
#                  axis.title.x = element_blank(), axis.title.y = element_blank())
# 
# p2 

p3 <- p1 + geom_hline(yintercept = -log10(SIG.LEVEL), linetype="longdash")
p3
p4 <- p3 + labs(x='Transcription Factor or Transcription Factor Pair', y='-log (pvalue)')

# p4 <- p3 + ggtitle(plot_title) + theme(plot.title = element_text(hjust = 0.5)) 
# 
# p4


# Save plots
#----------
# WIDTH <- 6
# HEIGHT <- 12
# fname.out <- paste(figdir, 'dotplots_logpvalues-model.props-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
# ggsave(filename = fname.out, p4, width = WIDTH, height = HEIGHT)

WIDTH <- 4
HEIGHT <- 8 #12
fname.out <- paste(figdir, 'stats.Chi-squared_test-', WIDTH, 'x', HEIGHT,'.pdf', sep = '') 
ggsave(filename = fname.out, p4, width = WIDTH, height = HEIGHT)

