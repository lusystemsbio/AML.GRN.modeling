remove(list = ls()) 

NUM.RACIPE.MODELS <- 2 #2000 #2 

# output directories
#-------------------- 
outdir <- './circuits.sim/'
dir.create(outdir)

# Load library and functions
#---------------------------
library(sRACIPE)
circuits <- readRDS(file = './circuits/circuits.uniq.rds')
length(circuits)
set.seed(100) 

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
rv = foreach(idx_name = names(circuits), .inorder = TRUE) %dopar% { 
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

