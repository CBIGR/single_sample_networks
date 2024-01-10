#!/usr/bin/Rscript

###################################################################################################################
#### Script to perform scaling on complete networks
###################################################################################################################

library(data.table)
library(dplyr)

base <- '/home/boris/Documents/PhD/single_sample/networks'
base_hpc <- '/user/gent/435/vsc43501/single_sample/network_outputs'
setwd(base)

scale_network <- function(networkfile, offset = FALSE){
  
  # networkfile <-'./LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet'
  # offset <- FALSE
  network <- fread(paste0(networkfile, '.csv'), data.table = FALSE, fill = TRUE, header = TRUE)
  if (offset == TRUE){
    network[,1] <- NULL
  }
  edges <- network[,c(1,2)]; weights <- network[,c(3:ncol(network))]
  
  # I will scale to values between -1 and 1 because CSN networks also have values between zero and 1 -> new approach now: use abs(max())
  
  #use x/max(of all samples) with absolute values
  max_weight <- max(abs(weights))
  weights <- weights/max_weight

  #apply(weights, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
  
  # maxs <- apply(weights, 2, max)
  # mins <- apply(weights, 2, min)
  # weights <- as.data.frame(scale(weights, center = (maxs + mins)/2, scale=(maxs-mins)/2))
  
  network_scaled <- cbind(edges, weights)
  write.table(network_scaled, file = paste0(networkfile, '_-1_1scaled.csv'), quote = FALSE, sep = '\t', row.names = FALSE)
}

scale_network('./LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet')
scale_network('./LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet')
scale_network('./SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet')
scale_network('./SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet')
scale_network('./CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet')
scale_network('./CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet')
scale_network('./ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet')
scale_network('./ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet')
scale_network('./SSPGI/SSPGI_brain_edges_2_75_HumanNet', offset = TRUE)
scale_network('./SSPGI/SSPGI_lung_edges_2_75_HumanNet', offset = TRUE)

# Also scale aggregate networks
scale_network('./aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet')
scale_network('./aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet')

