#!/usr/bin/Rscript


################################################################################################################
#### Top n connected genes as hubs + check enrichment for driver genes
################################################################################################################

library(plyr)
library(igraph)
library(data.table)
library(UpSetR)
library(reshape2)
library(ggplot2)
library(stringr)
library(ggvenn)
library(ggpubr)

################################################################################################################
#### Functions
################################################################################################################

convertToSymbol <- function(top200list){
  top200list <- str_replace(top200list, "\\s.*","")
  top200list <- str_replace(top200list, "\\(.*", "")
  return(top200list)
}

top200nodes <- function(network, number_of_drivers, offset=TRUE){
  # network <- paste0(base, "SSPGI/SSPGI_brain_edges_2_75_HumanNet_top_25000_edges.csv")
  # offset <- TRUE
  # number_of_drivers <- 5
  network <- fread(network, data.table=FALSE, header=TRUE, fill=TRUE)
  
  if (offset == TRUE){
    network[,1] <- NULL
  }
  n <- dim(network)[2]-2
  net <- network
  hub_list <- vector("list", length(colnames(network))-2)
  names(hub_list) <- colnames(network)[c(-1,-2)]
  
  for (i in 1:n){
    #i <- 2
    non_null_ind <- which(net[i+2] !=0)
    sample_netw <- net[non_null_ind, c(1,2,i+2)]
    connections <- as.data.frame(sort(table(c(sample_netw$reg, sample_netw$tar)), decreasing = TRUE))
    limit <- connections[number_of_drivers,]$Freq
    top200 <- connections[connections$Freq >= limit, ]
    top200 <- as.vector(top200$Var1)
    top200 <- convertToSymbol(top200)
    hub_list[[colnames(network)[i+2]]] <- top200
  }
  return(hub_list)
}

overlapWithDrivers <- function(top200list, cancerType, driverlist){
  
  # top200list <- SSN_brain
  # driverlist <- HGG_drivers
  # cancerType <- 'Glioblastoma'
  
  names(top200list) <- gsub('X', '', names(top200list))
  for (i in 1:length(names(top200list))){ # Convert top200 nodes to gene symbols
    top200list[[i]] <- convertToSymbol(top200list[[i]])
  }
  
  # Select samples belonging to the subtupe of interest, extract union of top200 nodes
  # sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,24)]
  # sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  samples <- sample_metadata[sample_metadata$Subtype == cancerType, ]
  top200list <- top200list[names(top200list) %in% samples$DepMap_ID]
  nodes <- unique(unlist(unname(top200list))) #Union of all driver genes found in samples of the subtype of interest
  
  
  if (cancerType == 'Medulloblastoma' | cancerType == 'Glioblastoma'){
    aggregate <- fread('./aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet.csv', fill = TRUE, header = TRUE, data.table = FALSE) #Should we use the complete aggregate network or the top 25k edges?
    aggregate$reg <- convertToSymbol(aggregate$reg)
    aggregate$tar <- convertToSymbol(aggregate$tar)
    HumanNet_genes <- unique(as.vector(as.matrix(aggregate[,c(1,2)])))
  } 
  
  # all_drivers <- fread('./Census_all-Jun_1_2022-10_08.tsv')[,1]
  # all_drivers <- all_drivers[all_drivers$`Gene Symbol` %in% HumanNet_genes, ] # Genes not present in the network cannot be identified as hubs, so leave those out
  gene_list <- list(All_genes = HumanNet_genes, Method_specific = nodes, Subtype_specific = driverlist$Symbol)
  overlap <- intersect(HumanNet_genes, intersect(nodes, driverlist$Symbol))
  print(overlap)
  # Create venn diagram
  plot <- ggvenn(gene_list, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"), show_percentage = FALSE, text_size = 4, set_name_size = 2.5)
  return(plot)
}

test_hyper <- function(method, driverlist, samples){
  
  names(method) <- gsub('X', '', names(method))
  method <- method[names(method) %in% samples$DepMap_ID]
  drivers_in_network <- intersect(driverlist$Symbol, HumanNet_genes)
  hubs <- unique(unlist(unname(method))) 
  # succes <- intersect(hubs, drivers_in_network)
  # possible_successes <- drivers_in_network
  # population_possible_sucesses <- length(HumanNet_genes) - length(drivers_in_network)
  #hubs
  
  # Calculate hypergeometric test p-value
  total_genes <- length(HumanNet_genes)  # Total number of genes in your dataset
  #total_genes <- length(all_drivers) # If you want to enrich against a background of known cancer driver genes instead of against the network
  driver_genes <- length(drivers_in_network)  # Number of known driver genes
  
  hub_gene_count <- length(hubs)  # Number of hub genes
  overlap_count <- length(intersect(hubs, drivers_in_network))  # Number of overlapping genes
  
  res <- list()
  
  # Calculate expected overlap by chance
  expected_overlap <- driver_genes * (hub_gene_count / total_genes)
  
  # Calculate fold change
  res$fold_change <- overlap_count / expected_overlap
  
  # Calculate hypergeometric test p-value
  res$test <- phyper(overlap_count - 1, driver_genes, total_genes - driver_genes, hub_gene_count, lower.tail = FALSE)
  
  #return(test)
  return(res)
}

################################################################################################################
#### Analysis
################################################################################################################


setwd('/home/boris/Documents/PhD/single_sample/networks')
base <- '/home/boris/Documents/PhD/single_sample/networks/'

# Now check enrichment for known drivers
MBL_drivers_intogen <- fread('./IntOGen-DriverGenes_MBL.tsv')[,1]
HGG_drivers_intogen <- fread('./IntOGen-DriverGenes_GBM.tsv')[,1]
MBL_drivers_census <- fread('./Census_MBL_19_06_2023.tsv')[,1]
HGG_drivers_census <- fread('./Census_GBM_19_06_2023.tsv')[,1]

MBL_drivers <- rbind(MBL_drivers_intogen, MBL_drivers_census, use.names=FALSE)
HGG_drivers <- rbind(HGG_drivers_intogen, HGG_drivers_census, use.names=FALSE)


all_drivers_intogen <- fread('./IntOGen-DriverGenes_all.tsv')[,1]
all_drivers_census <- fread('./Census_all_23_06_2023.tsv')[,1]
all_drivers <- rbind(all_drivers_intogen, all_drivers_census, use.names=FALSE)


aggregate <- fread('./aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet.csv', fill = TRUE, header = TRUE, data.table = FALSE) #Should we use the complete aggregate network or the top 25k edges?
aggregate$reg <- convertToSymbol(aggregate$reg)
aggregate$tar <- convertToSymbol(aggregate$tar)
HumanNet_genes <- unique(as.vector(as.matrix(aggregate[,c(1,2)])))
tmp1 <- intersect(HumanNet_genes, MBL_drivers$Symbol)
tmp2 <- intersect(HumanNet_genes, HGG_drivers$Symbol)
length(union(tmp1, tmp2))



# number_of_hubs <- seq(5,500,10)
#number_of_hubs <- seq(1,200,5)
number_of_hubs <- c(200)
p_values_HGG <- data.frame(matrix(nrow = 6, ncol=length(number_of_hubs)))
rownames(p_values_HGG) <- c('LIONESS_brain', 'SSN_brain', 'ssPCC_brain', 'CSN_brain', 'SSPGI_brain', 'SWEET_brain'); colnames(p_values_HGG) <- number_of_hubs
p_values_MBL <- data.frame(matrix(nrow = 6, ncol=length(number_of_hubs)))
rownames(p_values_MBL) <- c('LIONESS_brain', 'SSN_brain', 'ssPCC_brain', 'CSN_brain', 'SSPGI_brain', 'SWEET_brain'); colnames(p_values_MBL) <- number_of_hubs

fold_changes_HGG <- data.frame(matrix(nrow = 6, ncol=length(number_of_hubs)))
rownames(fold_changes_HGG) <- c('LIONESS_brain', 'SSN_brain', 'ssPCC_brain', 'CSN_brain', 'SSPGI_brain', 'SWEET_brain'); colnames(fold_changes_HGG) <- number_of_hubs
fold_changes_MBL <- data.frame(matrix(nrow = 6, ncol=length(number_of_hubs)))
rownames(fold_changes_MBL) <- c('LIONESS_brain', 'SSN_brain', 'ssPCC_brain', 'CSN_brain', 'SSPGI_brain', 'SWEET_brain'); colnames(fold_changes_MBL) <- number_of_hubs


aggregate <- fread('./aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet.csv', fill = TRUE, header = TRUE, data.table = FALSE) #Should we use the complete aggregate network or the top 25k edges?
aggregate$reg <- convertToSymbol(aggregate$reg)
aggregate$tar <- convertToSymbol(aggregate$tar)
HumanNet_genes <- unique(as.vector(as.matrix(aggregate[,c(1,2)])))


for (j in (1:length(number_of_hubs))){
  # j <- 1
  k <- number_of_hubs[j]
  # Find hub genes
  LIONESS_brain <- top200nodes(paste0(base, "LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), k, offset=TRUE)
  CSN_brain <- top200nodes(paste0(base, "CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv"), k, offset=FALSE)
  SSN_brain <- top200nodes(paste0(base,"SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), k, offset=TRUE)
  iENA_brain <- top200nodes(paste0(base,"ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), k, offset=TRUE)
  SSPGI_brain <- top200nodes(paste0(base,"SSPGI/SSPGI_brain_edges_2_75_HumanNet_top_25000_edges_symbol.csv"), k, offset=TRUE)
  SWEET_brain <- top200nodes(paste0(base,"SWEET/SWEET_brain_no_doubles_no_loops_HN_-1_1scaled_top_25000_edges.csv"), k, offset=FALSE)

  # # Read in sample metadata and group samples based on subtype
  sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,19,24)]
  sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  MBL_samples <- sample_metadata[sample_metadata$Subtype == 'Medulloblastoma', ]
  HGG_samples <- sample_metadata[sample_metadata$Subtype == 'Glioblastoma', ]
  
   #Providing cancer type as a variable doesn't work with paste0 so change this manually
  type <- 'Medulloblastoma'; print(type)
  SSN_MBL <- overlapWithDrivers(SSN_brain, type, MBL_drivers)
  LIONESS_MBL <- overlapWithDrivers(LIONESS_brain, type, MBL_drivers)
  SWEET_MBL <- overlapWithDrivers(SWEET_brain, type, MBL_drivers)
  ssPCC_MBL <- overlapWithDrivers(iENA_brain, type, MBL_drivers)
  CSN_MBL <- overlapWithDrivers(CSN_brain, type, MBL_drivers)
  SSPGI_MBL <- overlapWithDrivers(SSPGI_brain, type, MBL_drivers)
  
  
  # pdf(paste0('../results/top200_connected_nodes/brain/Venn/',type,'_top', k, 'nodes_driverOverlap.pdf'))
  # print(ggarrange(SSN_MBL, LIONESS_MBL, ssPCC_MBL, LIONESS_MBL, SSPGI_MBL, nrow=3, ncol=2, labels = c('SSN', 'LIONESS', 'iENA', 'CSN', 'SSPGI'), font.label = list(size = 8)))
  # dev.off()
  
  type <- 'Glioblastoma'; print(type)
  SSN_HGG <- overlapWithDrivers(SSN_brain, type, HGG_drivers)
  LIONESS_HGG <- overlapWithDrivers(LIONESS_brain, type, HGG_drivers)
  SWEET_HGG <- overlapWithDrivers(SWEET_brain, type, HGG_drivers)
  ssPCC_HGG <- overlapWithDrivers(iENA_brain, type, HGG_drivers)
  CSN_HGG <- overlapWithDrivers(CSN_brain, type, HGG_drivers)
  SSPGI_HGG <- overlapWithDrivers(SSPGI_brain, type, HGG_drivers)
  
  # pdf(paste0('../results/top200_connected_nodes/brain/Venn/',type,'_top', k, 'nodes_driverOverlap.pdf'))
  # print(ggarrange(SSN_HGG, LIONESS_HGG, ssPCC_HGG, CSN_HGG, SSPGI_HGG, nrow=3, ncol=2, labels = c('SSN', 'LIONESS', 'iENA', 'CSN', 'SSPGI'), font.label = list(size = 8)))
  # dev.off()
  
  p_values_MBL['LIONESS_brain', j] <- test_hyper(LIONESS_brain, MBL_drivers, MBL_samples)$test
  p_values_MBL['SSN_brain', j] <- test_hyper(SSN_brain,  MBL_drivers, MBL_samples)$test
  p_values_MBL['ssPCC_brain', j] <- test_hyper(iENA_brain, MBL_drivers, MBL_samples)$test
  p_values_MBL['CSN_brain', j] <- test_hyper(CSN_brain, MBL_drivers, MBL_samples)$test
  p_values_MBL['SSPGI_brain', j] <- test_hyper(SSPGI_brain, MBL_drivers, MBL_samples)$test
  p_values_MBL['SWEET_brain', j] <- test_hyper(SWEET_brain, MBL_drivers, MBL_samples)$test
  
  p_values_HGG['LIONESS_brain', j] <- test_hyper(LIONESS_brain, HGG_drivers, HGG_samples)$test
  p_values_HGG['SSN_brain', j] <- test_hyper(SSN_brain, HGG_drivers, HGG_samples)$test
  p_values_HGG['ssPCC_brain', j] <- test_hyper(iENA_brain, HGG_drivers, HGG_samples)$test
  p_values_HGG['CSN_brain', j] <- test_hyper(CSN_brain, HGG_drivers, HGG_samples)$test
  p_values_HGG['SSPGI_brain', j] <- test_hyper(SSPGI_brain, HGG_drivers, HGG_samples)$test
  p_values_HGG['SWEET_brain', j] <- test_hyper(SWEET_brain, HGG_drivers, HGG_samples)$test
  
  
  fold_changes_MBL['LIONESS_brain', j] <- test_hyper(LIONESS_brain, MBL_drivers, MBL_samples)$fold_change
  fold_changes_MBL['SSN_brain', j] <- test_hyper(SSN_brain,  MBL_drivers, MBL_samples)$fold_change
  fold_changes_MBL['ssPCC_brain', j] <- test_hyper(iENA_brain, MBL_drivers, MBL_samples)$fold_change
  fold_changes_MBL['CSN_brain', j] <- test_hyper(CSN_brain, MBL_drivers, MBL_samples)$fold_change
  fold_changes_MBL['SSPGI_brain', j] <- test_hyper(SSPGI_brain, MBL_drivers, MBL_samples)$fold_change
  fold_changes_MBL['SWEET_brain', j] <- test_hyper(SWEET_brain, MBL_drivers, MBL_samples)$fold_change
  
  fold_changes_HGG['LIONESS_brain', j] <- test_hyper(LIONESS_brain, HGG_drivers, HGG_samples)$fold_change
  fold_changes_HGG['SSN_brain', j] <- test_hyper(SSN_brain, HGG_drivers, HGG_samples)$fold_change
  fold_changes_HGG['ssPCC_brain', j] <- test_hyper(iENA_brain, HGG_drivers, HGG_samples)$fold_change
  fold_changes_HGG['CSN_brain', j] <- test_hyper(CSN_brain, HGG_drivers, HGG_samples)$fold_change
  fold_changes_HGG['SSPGI_brain', j] <- test_hyper(SSPGI_brain, HGG_drivers, HGG_samples)$fold_change
  fold_changes_HGG['SWEET_brain', j] <- test_hyper(SWEET_brain, HGG_drivers, HGG_samples)$fold_change
}

save.image('../results/Rebuttal/top200_connected_nodes/results_top200_brain_Intogen_drivers_withSWEET.RData')
#load('../results/top200_connected_nodes/brain/results5_500.RData')

names(SSN_brain) <- gsub('X','', names(SSN_brain))
names(LIONESS_brain) <- gsub('X','', names(LIONESS_brain))
names(iENA_brain) <- gsub('X','', names(iENA_brain))
names(CSN_brain) <- gsub('X','', names(CSN_brain))
names(SSPGI_brain) <- gsub('X','', names(SSPGI_brain))
names(SWEET_brain) <- gsub('X','', names(SWEET_brain))

# Define your sample names
sample_names <- names(LIONESS_brain)

# Create an empty data frame to store the overlap results
overlap_results <- list()

# Loop through each sample and calculate the overlap
for (sample in sample_names) {
  
  #sample <- '0010'
  
  # Get the sets of unique values for each method
  lioness_set <- unlist(LIONESS_brain[[sample]])
  ssn_set <- unlist(SSN_brain[[sample]])
  sweet_set <- unlist(SWEET_brain[[sample]])
  iena_set <- unlist(iENA_brain[[sample]])
  csn_set <- unlist(CSN_brain[[sample]])
  sspgi_set <- unlist(SSPGI_brain[[sample]])
  
  # Calculate the overlap across all 6 lists
  #overlap_all_methods <- length(Reduce(intersect, list(lioness_set, ssn_set, iena_set, csn_set, sspgi_set)))
  overlap_all_methods <- length(Reduce(intersect, list(lioness_set, ssn_set, sweet_set, iena_set, csn_set, sspgi_set)))
  overlap_results[[sample]] <- overlap_all_methods
}
# Display the results for each sample
for (sample in sample_names) {
  cat("Sample:", sample, "\n")
  cat("Overlapping Features:", overlap_results[[sample]], "\n")
}
# Display the overlap results
#print(overlap_results)

# Calculate the average overlap for each method
average_overlap <- sapply(overlap_results, mean)
mean(unlist(average_overlap))

LIONESS_intersection <- Reduce(intersect, unname(LIONESS_brain)) # zero genes overlapping across all samples
SSN_intersection <- Reduce(intersect, SSN_brain) # zero genes
CSN_intersection <- Reduce(intersect, CSN_brain) # 76 genes (so really a lot if you think about the fact that we use top200 as hubs)
iENA_intersection <- Reduce(intersect, iENA_brain) # zero genes
SSPGI_intersection <- Reduce(intersect, SSPGI_brain) # 21 genes
