
remove(list = ls())

source('./functions.tfBalance.R')

inputDir <- './data/'

fname.data <- paste(inputDir, 'stats.upTF.csv', sep = '')
stats.upTF.df <- read.csv(file = fname.data, header = T) 
YLIMIT <- c(0, 44)
XLIMIT <- c(1, 170)

figdir <- './figs/'
dir.create(figdir)

method <- "Netact"
stats.df <- stats.upTF.df[stats.upTF.df$method==method, ]
figname <- paste(figdir, 'stats.', method, '.pdf', sep = '') 
plotStat(stats.df=stats.df, colors, figname=figname, 
         ylimit=YLIMIT, xlimit=XLIMIT)

method <- "MARINa"
stats.df <- stats.upTF.df[stats.upTF.df$method==method, ]
figname <- paste(figdir, 'stats.', method, '.pdf', sep = '') 
plotStat(stats.df=stats.df, colors, figname=figname, 
         ylimit=YLIMIT, xlimit=XLIMIT)

method <- "RI"
stats.df <- stats.upTF.df[stats.upTF.df$method==method, ]
figname <- paste(figdir, 'stats.', method, '.pdf', sep = '') 
plotStat(stats.df=stats.df, colors, figname=figname, 
         ylimit=YLIMIT, xlimit=XLIMIT)

method <- "COMB"
stats.df <- stats.upTF.df[stats.upTF.df$method==method, ]
figname <- paste(figdir, 'stats.', method, '.pdf', sep = '') 
plotStat.comb(stats.df=stats.df, colors, figname=figname, 
              ylimit=YLIMIT*3, xlimit=XLIMIT)


colors <- c('black', 'magenta','blue', 'red') 
names(colors) <- c("TFs", "ACTUAL", "CTRL", "TRTMT") 
ylimit=c(0, 150)
plot(stats.df$no_tfs*3, col=colors["TFs"], pch = 19,  
     xlim = c(1, 110), ylim = ylimit, xaxt = 'none',
     xlab = 'feature ratio', ylab = 'Number of TFs')
points(stats.df$no_tfs_actual, col=colors["ACTUAL"])

plot(stats.df$no_tfs*3-stats.df$no_tfs_actual)

