find_TFs.R
   inputs: 
      (1) expression data: raw_brainarray.sele.txt 
      (2) lib.aml.idh.R
   Operations: 
      (1) creates Expression Set object and writes it out to the 
          file eset.brain_array.rda 
      (2) find DE genes for each comparison and writes it to the file 
          de.results.rda 
      (3) find ranked TFs and their statistics for each comparison group
   outputs: 
      (1) eset.brain_array.rda
      (2) de.results.rda
      (3) tfs.CTRL.IDH.csv
          contains TFs from TF_GSEA for group CTRL-IDH
         
order.TFs
   inputs: 
      (1) TFs: tfs.CTRL.IDH.csv
   Operations: for each CTRL-TRMT group -  
      (1) order TFs whose Qvalues are zero by their Z scores
      (2) create a new column zq where 
          Z scores are found for the TFs whose qvalues are zero
          Q values around for the subsequent TFs
   outputs: 
      ordered TFs with the new column zq: tfs.CTR.IDH.ordered.csv

select.TFs.R
   inputs: 
      (1) read ordered TFs files: tfs.ordered.rds
   Operations: for each comparison group -  
      (1) Select balanced number of TFs from CTRL and TRTMT
      (2) Filter TOP TFs based on Pvalue score from each group
      (2) merge selected TFs
      (3) save merged TFs
   outputs:
      (1) TOP TFs obtained after merging: tfs.TOP.phase.31.txt

comb.TOP.TFs.R
   inputs: 
      (1) tfs.TOP.phase.31.txt 
      (2) ./from.xiaowen/masterTFs.MARINa.txt
   operations: combined both TF sets 
   outputs: tfs.TOP.txt

cal.ACT_EXP.R
   inuts: 
      (1) TOP tfs: tfs.TOP.txt
   operations:
      (1) calculate activities by invoking TF_Activity for the TOP tfs 
      (2) plot activities vs expressions heatmaps by invoking Combine_heatmap
   outputs: 
      (1) combined heatmaps

cal.ACT_EXP.diagnose.R
   inputs:  
      (1) tfs.TOP.txt 
      (2) amlnet.tpo - 26x72 - network used in phase.13 
   operations: 
      (1) 
       
   outputs: 
    figs.diagnose

