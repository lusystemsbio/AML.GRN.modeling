remove(list = ls()) 

NUM.RACIPE.MODELS <- 100000 #2 
NO_TOP_CIRCUITS <- 11

# output directories
#-------------------- 
outdir <- './circuits.sim/'
dir.create(outdir)

# Load library and functions
#---------------------------
library(sRACIPE)

# Load circuits
#-------------- 
circuits <- readRDS(file = '../../networks/circuits/circuits.uniq.rds')
length(circuits)
set.seed(100) 

# Load circuit metrics
#---------------------------
#circuit_metrics.sim <- read.csv(file = '../networks/results/summary.circuits.sim.sorted.csv', row.names = 1)
circuit_metrics.sim <- read.csv(file = '../rankings.mean_acc/results/summary.circuits.sim.sortedByAcc_flex.csv', row.names = 1)
topCircuits <- rownames(circuit_metrics.sim)[1:NO_TOP_CIRCUITS]

# Allocate resources for parallel computing
#-------------------------------------------#
library(doParallel)
library(foreach)
NO_AVAILABLE_CORES <- detectCores() 
print(NO_AVAILABLE_CORES)
cl <- makeCluster(NO_AVAILABLE_CORES)
class(cl)
registerDoParallel(cl)
getDoParWorkers()


# Execute RACIPE simulations  - parallel
#----------------------------------------#
rv = foreach(idx_name = topCircuits, .inorder = TRUE) %dopar% { 
      library(sRACIPE)
      circuitSimulated <- sracipeSimulate(circuit = circuits[[idx_name]][,c("Source", "Target", "Type")], 
                                  plotToFile = TRUE,
                                  numModels = NUM.RACIPE.MODELS, # Change 
                                  plots = FALSE,
                                  stepper = "RK4",
                                  integrateStepSize = 0.05) 
      fname.sim.circuit <- paste(outdir, './circuit_simulated_', idx_name ,'.rds', sep = '') 
      saveRDS(circuitSimulated, file = fname.sim.circuit)
}
 
#  Release resources for parallel computation 
#---------------------------------------------#
stopCluster(cl)

