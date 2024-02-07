
# Create mapping between top circuit ids and their names (given as circuit_1, ...)
create_circuit_id_2_name_mapping <- function(top.circuit_id.list){
   top.circuits <- unique(as.character(unlist(top.circuit_id.list))) 
   names(top.circuits) <- paste('circuit_', 1:length(top.circuits), sep = '')
   
   mapping_circuit_ids <- as.data.frame(cbind(as.character(top.circuits), names(top.circuits))) 
   colnames(mapping_circuit_ids) <- c('circuit_id', 'circuit_name')
   rownames(mapping_circuit_ids) <- names(top.circuits)  
   return(mapping_circuit_ids)
}



get_barPlot <- function(mean_sdTOPmodels){ 
   library(ggplot2)
   p1 <- ggplot(mean_sdTOPmodels) + 
      geom_bar( aes(x=circuit_id, y=mean), 
                stat="identity", fill="forestgreen", alpha=0.5) +
      geom_errorbar(aes(x=circuit_id, 
                        ymin=mean-sd, ymax=mean+sd), width=0.4, 
                    colour="orange", alpha=0.9, size=1.5) #+
      #ggtitle("mean of SDs")
   
   #p1
   
   p2 <- p1 + theme(axis.text.x = element_text(angle = 90), 
                    axis.title.x = element_blank(), 
                    axis.title.y = element_blank())
   #p2   
   return(p2)
}
