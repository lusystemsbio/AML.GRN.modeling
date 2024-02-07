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
simulatedCircuits <- list()
for(idx_name in names(circuits)){
      simulatedCircuits[[idx_name]] <- assay(sracipeSimulate(circuit = circuits[[idx_name]][,c("Source", "Target", "Type")], 
                                                       plotToFile = TRUE,
                                                       numModels = NUM.RACIPE.MODELS, #20,  # Change 
                                                       plots = FALSE, 
                                                       stepper = "RK4", 
                                                       integrateStepSize = 0.05)
                                             )
      #simulatedCircuits[[idx_name]] <- simulatedCircuit
   #break()
}

fname.sim.circuits <- paste0(outdir, './circuits.simulated.rds')
saveRDS(simulatedCircuits, file = fname.sim.circuits)


