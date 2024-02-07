#!/usr/bin/env Rscript
remove(list = ls())

circuit_metrics.sim <- read.csv(file = './results/summary.circuits.sim.csv', row.names = 1)

flexibility.df <- read.csv(file = './results/circuit.flexibility.csv')
rownames(flexibility.df) <- flexibility.df$circuit.idx

# scale flexibility 
# accuracy.stat <- mean(circuit_metrics.sim$Accuracy)
# flex.stat <- mean(flexibility.df$flexibility)
# acc_flex_ratio <- accuracy.stat/flex.stat
# acc_flex_ratio # 17.14
# 
# accuracy.stat <- median(circuit_metrics.sim$Accuracy)
# flex.stat <- median(flexibility.df$flexibility)
# acc_flex_ratio <- accuracy.stat/flex.stat
# acc_flex_ratio


# tmp.df <- flexibility.df 
# flexibility.df["flexibility.scaled"] <- flexibility.df$flexibility*acc_flex_ratio 
# 
# flexibility.df["flexibility.scaled"] <- (flexibility.df["flexibility.scaled"]-min(flexibility.df["flexibility.scaled"]))/(max(flexibility.df["flexibility.scaled"])-min(flexibility.df["flexibility.scaled"]))
# plot(flexibility.df$flexibility.scaled)

flexibility.ordered <- flexibility.df[rownames(circuit_metrics.sim), ]

sum(rownames(circuit_metrics.sim) == rownames(flexibility.ordered))

circuit_metrics.sim.new <- cbind(circuit_metrics.sim, flexibility.ordered$flexibility)
colnames(circuit_metrics.sim.new) <- c(colnames(circuit_metrics.sim), 'flexibility')

circuit_metrics.sim.new <- cbind(circuit_metrics.sim, flexibility.ordered$flexibility)
colnames(circuit_metrics.sim.new) <- c(colnames(circuit_metrics.sim), 'flexibility')

outdir <- './results/'
dir.create(outdir)
write.csv(circuit_metrics.sim.new, file = paste0(outdir, "./summary.circuits.sim.flex.csv"),
          row.names = T, quote = F)
