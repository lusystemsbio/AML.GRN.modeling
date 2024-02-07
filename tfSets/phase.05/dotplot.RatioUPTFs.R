
remove(list = ls())

source('./functions.tfBalance.R')
inputDir <- './data/'

fname.data <- paste(inputDir, 'stats.upTF.csv', sep = '')
stats.upTF.df <- read.csv(file = fname.data, header = T) 

min(stats.upTF.df$TRTMT_CTRL)
max(stats.upTF.df$TRTMT_CTRL)

XLIMIT <- c(1, 170)
YLIMIT <- c(0, 1)

figdir <- './figs/'
dir.create(figdir)

figname <- paste(figdir, 'ratio-upTFs.', '.pdf', sep = '')
WIDTH <- 10
HEIGHT <- 8 
pdf(file = figname, width=WIDTH, height=HEIGHT, paper='special')

par(mfrow=c(4,1))
par(mar=c(1.0, 3.1, 1.1, 1.1), # bottom, left, top, right
    mgp=c(3, 1, 0), las=0)

method <- "Netact"
ratio_upTFs <- stats.upTF.df[stats.upTF.df$method==method, 'TRTMT_CTRL']
plot(stats.upTF.df$TRTMT_CTRL, pch = 19, type = 'n',
     xlim = XLIMIT, ylim = YLIMIT, xaxt='n')
      
points(ratio_upTFs)

method <- "MARINa"
ratio_upTFs <- stats.upTF.df[stats.upTF.df$method==method, 'TRTMT_CTRL']
plot(stats.upTF.df$TRTMT_CTRL, pch = 19, type = 'n',
     xlim = XLIMIT, ylim = YLIMIT, xaxt='n') 
points(ratio_upTFs)
 
method <- "RI"
ratio_upTFs <- stats.upTF.df[stats.upTF.df$method==method, 'TRTMT_CTRL']
plot(stats.upTF.df$TRTMT_CTRL, pch = 19, type = 'n',
     xlim = XLIMIT, ylim = YLIMIT, xaxt='n') 
points(ratio_upTFs)

method <- "COMB"
ratio_upTFs <- stats.upTF.df[stats.upTF.df$method==method, 'TRTMT_CTRL']
plot(stats.upTF.df$TRTMT_CTRL, pch = 19, type = 'n',
     xlim = XLIMIT, ylim = YLIMIT, xaxt='n')
points(ratio_upTFs) 
dev.off()

