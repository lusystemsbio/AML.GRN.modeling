#/usr/bin/env bash 

# obtain DEGs between AML and CTRL data: 
#R CMD BATCH obtain.DEGs.R

# obtain DEGs that are either circuit TFs or their targets:
#R CMD BATCH obtain.DEGs_in_circuit.R

# obtain significan pathways using TF_GSEA function from NetAct 
# OR
#R CMD BATCH pathway.analysis.GSEA.R
# obtain significant pathways using enrichr 

# Fisher's exact test between each pathway pair: 
# calculate relation between pathways - we DO NOT need this step:
#R CMD BATCH cal.overlap_btn_pathways.R 


# map TFs to pathways:
# (1) for the qualified TFs, assign significant pathways
# (2) for each tf, assign top significant pathway when available,
#     otherwise assign NA
# (3) create a mapping between pathway name to pathway id
R CMD BATCH map_TFs_2_pathways.R

# extract and save circuit topology for the top circuit
R CMD BATCH xtract.topology.R
