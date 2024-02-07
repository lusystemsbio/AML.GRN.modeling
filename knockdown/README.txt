current folder (./), ./data, ./figs: 
 activities to combine nnet predictions for skd and dkd 

tpo/
 circuit topology file

racipe.wt/
 Untreated RACIPE models

skd/
 activities performed on single knockdwon(SKD) racipe models

dkd/
 activities performed on double knockdwon(SKD) racipe models

cluster_props.comb/
 combine the cluster proportions from skd and dkd 
 plot figures for the combined cluster proportions

nnet.performance/
 The proportions of models from single and double knockdown are 
 combined at different probability difference threshold level 
 and plotted.

nnet.robustness/
   the robustness (sensitivity vs specificity) of nnet prediction model 
   is evaluated. 
   
chiSquared.test/
   Chi-squared test is performed on the combined (single and double knockdown) 
   proportion of models.  

slides/

wip/

# Depracated
-------------
comb.hS/
 activities to combine heatmapSimilarity predictions for skd and dkd 

