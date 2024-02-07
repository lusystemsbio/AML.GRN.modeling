remove(list = ls()) 

NUM.RACIPE.MODELS <- 2 #2000 #2 

# output directories
#-------------------- 
outdir <- './circuits/'


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
simulatedCircuits <- foreach(idx_name = names(circuits), .inorder = TRUE) %dopar% { 
   library(sRACIPE)
   assay(sracipeSimulate(circuit = circuits[[idx_name]][,c("Source", "Target", "Type")], 
                   plotToFile = TRUE,
                   numModels = NUM.RACIPE.MODELS, # Change 
                   plots = FALSE,
                   stepper = "RK4",
                   integrateStepSize = 0.05)
         )
}
names(simulatedCircuits) <- names(circuits)

#  Release resources for parallel computation 
#---------------------------------------------#
stopCluster(cl)


# Save 
#--------#
fname.sim.circuits <- paste0(outdir, './circuits.simulated.rds')
saveRDS(simulatedCircuits, file = fname.sim.circuits)

