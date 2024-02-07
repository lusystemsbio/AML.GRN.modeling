




circuit.nodes.Netact <- coreTFs.Netact[which(coreTFs.Netact %in% circuit.nodes)] 

sigvalsNetact <- c(min(as.numeric(labels(circuit.nodes.Netact))), 
                   max(as.numeric(labels(circuit.nodes.Netact))))
sigvalsNetact


# obtain signifiance values of the TFs
coreTFs.Method <- coreTFs.Netact







