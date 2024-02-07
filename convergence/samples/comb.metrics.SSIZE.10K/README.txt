Step 1. Calculate metrics across samples S1 ~ S10 
    metrics: Accuracy, AvgDist, and Flexibility 
 cal.metrics.sim_circuits.R

Step 2. Rank the circuits: 3-way sorting using Accuracy and Flexibility 
 rankCircuits.R

Step 3. Calculate mean and sd of combined scores (ranks) across samples for each circuit: 
 cal.mean.sd.R
