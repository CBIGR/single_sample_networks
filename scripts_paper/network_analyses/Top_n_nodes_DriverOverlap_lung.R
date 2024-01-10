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
  
  # network <- paste0(base, "SSPGI/SSPGI_lung_edges_2_75_HumanNet_top_25000_edges_symbol.csv")
  # offset <- TRUE
  # number_of_drivers <- 200
  
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
  
  # top200list <- LIONESS_lung
  # driverlist <- Squamous_drivers
  # cancerType <- 'NSCLC_squamous'
  
  names(top200list) <- gsub('X', '', names(top200list))
  for (i in 1:length(names(top200list))){ # Convert top200 nodes to gene symbols
    top200list[[i]] <- convertToSymbol(top200list[[i]])
  }
  
  # Select samples belonging to the subtupe of interest, extract union of top200 nodes
  # sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,24)]
  # sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  samples <- sample_metadata[sample_metadata$lineage_sub_subtype == cancerType, ]
  top200list <- top200list[names(top200list) %in% samples$DepMap_ID]
  nodes <- unique(unlist(unname(top200list))) #Union of all driver genes found in samples of the subtype of interest
  
  
  if (cancerType == 'SCLC' | cancerType == 'NSCLC_adenocarcinoma' | cancerType == 'NSCLC_squamous'){
    aggregate <- fread('./aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv', fill = TRUE, header = TRUE, data.table = FALSE) #Should we use the complete aggregate network or the top 25k edges?
    aggregate$reg <- convertToSymbol(aggregate$reg)
    aggregate$tar <- convertToSymbol(aggregate$tar)
    HumanNet_genes <- unique(as.vector(as.matrix(aggregate[,c(1,2)])))
  } 
  
  
  # all_drivers <- fread('./Census_all-Jun_1_2022-10_08.tsv')[,1]
  # all_drivers <- all_drivers[all_drivers$`Gene Symbol` %in% HumanNet_genes, ] # Genes not present in the network cannot be identified as hubs, so leave those out
  gene_list <- list(All_genes = HumanNet_genes, Method_specific = nodes, Subtype_specific = driverlist$Symbol)
  
  # Create venn diagram
  overlap <- intersect(HumanNet_genes, intersect(nodes, driverlist$Symbol))
  print(overlap)
  plot <- ggvenn(gene_list, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"), show_percentage = FALSE, text_size = 4, set_name_size = 2.5)
  return(plot)
}

test_hyper <- function(method, driverlist, samples){
 
  # methhod <- SSN_lung
  # samples <- SCLC_samples
  # driverlist <- SCLC_drivers
  
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

aggregate <- fread('./aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv', fill = TRUE, header = TRUE, data.table = FALSE) #Should we use the complete aggregate network or the top 25k edges?
aggregate$reg <- convertToSymbol(aggregate$reg)
aggregate$tar <- convertToSymbol(aggregate$tar)

# Now check enrichment for known drivers
SCLC_drivers_intogen <- fread('./IntOGen-DriverGenes_SCLC.tsv')[,1]
#NSCLC_drivers_intogen <- fread('./IntOGen-DriverGenes_NSCLC.tsv')[,1]
SCLC_drivers_census <- fread('./Census_SCLC_19_06_2023.tsv')[,1]
#NSCLC_drivers_census <- fread('./Census_NSCLC_19_06_2023.tsv')[,1]
Adeno_drivers_intogen <- fread('./IntOGen-DriverGenes_Lung-Adenocarcinoma.tsv')[,1]
Adeno_drivers_census <- fread('./Census_Adenocarcinoma_16_06_2023.tsv')[,1]
Squamous_drivers_intogen <- fread('./IntOGen-DriverGenes_Lung-Squamous.tsv')[,1]
NSCLC_drivers_intogen <- fread('./IntOGen-DriverGenes_NSCLC.tsv')[,1]
NSCLC_drivers_census <- fread('./Census_NSCLC_19_06_2023.tsv')[,1]



SCLC_drivers <- rbind(SCLC_drivers_intogen, SCLC_drivers_census, use.names=FALSE)
NSCLC_drivers <- rbind(NSCLC_drivers_intogen, NSCLC_drivers_census, use.names=FALSE)
Adeno_drivers <- rbind(Adeno_drivers_census, Adeno_drivers_intogen, use.names = FALSE); colnames(Adeno_drivers) <- c('Symbol')
Squamous_drivers <- Squamous_drivers_intogen
#NSCLC_drivers <- rbind(NSCLC_drivers_intogen, NSCLC_drivers_census, use.names=FALSE)

all_drivers_intogen <- fread('./IntOGen-DriverGenes_all.tsv')[,1]
all_drivers_census <- fread('./Census_all_23_06_2023.tsv')[,1]
all_drivers <- rbind(all_drivers_intogen, all_drivers_census, use.names=FALSE)


# HumanNet_genes <- unique(as.vector(as.matrix(aggregate[,c(1,2)])))
# tmp1 <- intersect(HumanNet_genes, Adeno_drivers$Symbol)
# tmp2 <- intersect(HumanNet_genes, SCLC_drivers$Symbol)
# tmp3 <- intersect(HumanNet_genes, Squamous_drivers$Symbol)
#length(union(tmp1, tmp2, tmp3))


# number_of_hubs <- seq(5,500,10)
#number_of_hubs <- seq(1,200,5)
number_of_hubs <- c(200)
p_values_SCLC <- data.frame(matrix(nrow = 6, ncol=length(number_of_hubs)))
rownames(p_values_SCLC) <- c('LIONESS_lung', 'SSN_lung', 'ssPCC_lung', 'CSN_lung', 'SSPGI_lung', 'SWEET_lung'); colnames(p_values_SCLC) <- number_of_hubs
p_values_Adeno <- data.frame(matrix(nrow = 6, ncol=length(number_of_hubs)))
rownames(p_values_Adeno) <- c('LIONESS_lung', 'SSN_lung', 'ssPCC_lung', 'CSN_lung', 'SSPGI_lung', 'SWEET_lung'); colnames(p_values_Adeno) <- number_of_hubs
p_values_Squamous <- data.frame(matrix(nrow = 6, ncol=length(number_of_hubs)))
rownames(p_values_Squamous) <- c('LIONESS_lung', 'SSN_lung', 'ssPCC_lung', 'CSN_lung', 'SSPGI_lung', 'SWEET_lung'); colnames(p_values_Squamous) <- number_of_hubs
p_values_NSCLC <- data.frame(matrix(nrow = 6, ncol=length(number_of_hubs)))
rownames(p_values_NSCLC) <- c('LIONESS_lung', 'SSN_lung', 'ssPCC_lung', 'CSN_lung', 'SSPGI_lung', 'SWEET_lung'); colnames(p_values_NSCLC) <- number_of_hubs

fold_changes_SCLC <- data.frame(matrix(nrow = 6, ncol=length(number_of_hubs)))
rownames(fold_changes_SCLC) <- c('LIONESS_lung', 'SSN_lung', 'ssPCC_lung', 'CSN_lung', 'SSPGI_lung', 'SWEET_lung'); colnames(fold_changes_SCLC) <- number_of_hubs
fold_changes_Adeno <- data.frame(matrix(nrow = 6, ncol=length(number_of_hubs)))
rownames(fold_changes_Adeno) <- c('LIONESS_lung', 'SSN_lung', 'ssPCC_lung', 'CSN_lung', 'SSPGI_lung', 'SWEET_lung'); colnames(fold_changes_Adeno) <- number_of_hubs
fold_changes_Squamous <- data.frame(matrix(nrow = 6, ncol=length(number_of_hubs)))
rownames(fold_changes_Squamous) <-c('LIONESS_lung', 'SSN_lung', 'ssPCC_lung', 'CSN_lung', 'SSPGI_lung', 'SWEET_lung'); colnames(fold_changes_Squamous) <- number_of_hubs
fold_changes_NSCLC <- data.frame(matrix(nrow = 6, ncol=length(number_of_hubs)))
rownames(fold_changes_NSCLC) <- c('LIONESS_lung', 'SSN_lung', 'ssPCC_lung', 'CSN_lung', 'SSPGI_lung', 'SWEET_lung'); colnames(fold_changes_NSCLC) <- number_of_hubs


aggregate <- fread('./aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv', fill = TRUE, header = TRUE, data.table = FALSE) #Should we use the complete aggregate network or the top 25k edges?
aggregate$reg <- convertToSymbol(aggregate$reg)
aggregate$tar <- convertToSymbol(aggregate$tar)
HumanNet_genes <- unique(as.vector(as.matrix(aggregate[,c(1,2)])))


# number_of_hubs <- c(5, 200)
for (j in (1:length(number_of_hubs))){
  k <- 200
  
  # k <- number_of_hubs[j]
  # Find hub genes
  LIONESS_lung <- top200nodes(paste0(base, "LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv"),k, offset=TRUE)
  CSN_lung <- top200nodes(paste0(base, "CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv"), k, offset=FALSE)
  SSN_lung <- top200nodes(paste0(base,"SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), k, offset=TRUE)
  iENA_lung <- top200nodes(paste0(base,"ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), k, offset=TRUE)
  SSPGI_lung <- top200nodes(paste0(base,"SSPGI/SSPGI_lung_edges_2_75_HumanNet_top_25000_edges_symbol.csv"), k, offset=TRUE)
  SWEET_lung <- top200nodes(paste0(base,"SWEET/SWEET_lung_no_doubles_no_loops_HN_-1_1scaled_top_25000_edges.csv"), k, offset=TRUE)
  # save.image(paste0('../results/top200_connected_nodes/lung/top',k,'nodes.RData'))
  
  # Read in sample metadata and group samples based on subtype
  sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,24,25)]
  sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  # SCLC samples have no annotation in sub_subtype, just put the subtype there
  sample_metadata$lineage_sub_subtype[sample_metadata$lineage_sub_subtype == ""] <- sample_metadata$lineage_subtype[sample_metadata$lineage_sub_subtype == ""]

  
  SCLC_samples <- sample_metadata[sample_metadata$lineage_sub_subtype == 'SCLC', ]
  Adeno_samples <- sample_metadata[sample_metadata$lineage_sub_subtype == 'NSCLC_adenocarcinoma']
  Squamous_samples <- sample_metadata[sample_metadata$lineage_sub_subtype == 'NSCLC_squamous']
  LargeCell_samples <- sample_metadata[sample_metadata$lineage_sub_subtype == 'NSCLC_large_cell']
  NSCLC_samples <- sample_metadata[sample_metadata$lineage_subtype == 'NSCLC']
  
  # Once for NSCLC and once for SCLC
  type <- 'SCLC'
  print('SCLC overlap with known drivers')
  
  SSN_SCLC <- overlapWithDrivers(SSN_lung, type, SCLC_drivers)
  LIONESS_SCLC <- overlapWithDrivers(LIONESS_lung, type, SCLC_drivers)
  SWEET_SCLC <- overlapWithDrivers(SWEET_lung, type, SCLC_drivers)
  ssPCC_SCLC <- overlapWithDrivers(iENA_lung, type, SCLC_drivers)
  CSN_SCLC <- overlapWithDrivers(CSN_lung, type, SCLC_drivers)
  SSPGI_SCLC <- overlapWithDrivers(SSPGI_lung, type, SCLC_drivers)
  
  
  type <- 'NSCLC'
  print('NSCLC overlap with known drivers')
  
  SSN_NSCLC <- overlapWithDrivers(SSN_lung, type, NSCLC_drivers)
  LIONESS_NSCLC <- overlapWithDrivers(LIONESS_lung, type, NSCLC_drivers)
  SWEET_NSCLC <- overlapWithDrivers(SWEET_lung, type, NSCLC_drivers)
  ssPCC_NSCLC <- overlapWithDrive#!/usr/bin/Rscript
rs(iENA_lung, type, NSCLC_drivers)
  CSN_NSCLC <- overlapWithDrivers(CSN_lung, type, NSCLC_drivers)
  SSPGI_NSCLC <- overlapWithDrivers(SSPGI_lung, type, NSCLC_drivers)
  
  
  # pdf(paste0('../results/top200_connected_nodes/lung/Venn/',type,'_top', k, 'nodes_driverOverlap.pdf'))
  # print(ggarrange(SSN_SCLC, LIONESS_SCLC, ssPCC_SCLC, CSN_SCLC, SSPGI_SCLC, nrow=3, ncol=2, labels = c('SSN', 'LIONESS', 'iENA', 'CSN', 'SSPGI'), font.label = list(size = 8)))
  # dev.off()
  
  # Once for NSCLC and once for SCLC
  type <- 'NSCLC_adenocarcinoma'
  print('Adenocarcinoma')
  SSN_Adeno <- overlapWithDrivers(SSN_lung, type, Adeno_drivers)
  LIONESS_Adeno <- overlapWithDrivers(LIONESS_lung, type, Adeno_drivers)
  SWEET_Adeno <- overlapWithDrivers(SWEET_lung, type, Adeno_drivers)
  ssPCC_Adeno <- overlapWithDrivers(iENA_lung, type, Adeno_drivers)
  CSN_Adeno <- overlapWithDrivers(CSN_lung, type, Adeno_drivers)
  SSPGI_Adeno <- overlapWithDrivers(SSPGI_lung, type, Adeno_drivers)
  
  
  # Once for NSCLC and once for SCLC
  type <- 'NSCLC_squamous'
  print('Squamous')
  SSN_Squamous <- overlapWithDrivers(SSN_lung, type, Squamous_drivers)
  LIONESS_Squamous <- overlapWithDrivers(LIONESS_lung, type, Squamous_drivers)
  SWEET_Squamous <- overlapWithDrivers(SWEET_lung, type, Squamous_drivers)
  ssPCC_Squamous <- overlapWithDrivers(iENA_lung, type, Squamous_drivers)
  CSN_Squamous <- overlapWithDrivers(CSN_lung, type, Squamous_drivers)
  SSPGI_Squamous <- overlapWithDrivers(SSPGI_lung, type, Squamous_drivers)
  
  
  # pdf(paste0('../results/top200_connected_nodes/lung/Venn/',type,'_top', k, 'nodes_driverOverlap.pdf'))
  # print(ggarrange(SSN_NSCLC, LIONESS_NSCLC, ssPCC_NSCLC, CSN_NSCLC, SSPGI_NSCLC, nrow=3, ncol=2, labels = c('SSN', 'LIONESS', 'iENA', 'CSN', 'SSPGI'), font.label = list(size = 8)))
  # dev.off()
  
  # Now test enrichment under a hypergeometric distribution for both SCLC and NSCLC
  p_values_SCLC['LIONESS_lung', j] <- test_hyper(LIONESS_lung, SCLC_drivers, SCLC_samples)$test
  p_values_SCLC['SSN_lung', j] <- test_hyper(SSN_lung, SCLC_drivers, SCLC_samples)$test
  p_values_SCLC['ssPCC_lung', j] <- test_hyper(iENA_lung, SCLC_drivers, SCLC_samples)$test
  p_values_SCLC['CSN_lung', j] <- test_hyper(CSN_lung, SCLC_drivers, SCLC_samples)$test
  p_values_SCLC['SSPGI_lung', j] <- test_hyper(SSPGI_lung, SCLC_drivers, SCLC_samples)$test
  p_values_SCLC['SWEET_lung', j] <- test_hyper(SWEET_lung, SCLC_drivers, SCLC_samples)$test
  
  p_values_Adeno['LIONESS_lung', j] <- test_hyper(LIONESS_lung, Adeno_drivers, Adeno_samples)$test
  p_values_Adeno['SSN_lung', j] <- test_hyper(SSN_lung, Adeno_drivers, Adeno_samples)$test
  p_values_Adeno['ssPCC_lung', j] <- test_hyper(iENA_lung, Adeno_drivers, Adeno_samples)$test
  p_values_Adeno['CSN_lung', j] <- test_hyper(CSN_lung, Adeno_drivers, Adeno_samples)$test
  p_values_Adeno['SSPGI_lung', j] <- test_hyper(SSPGI_lung, Adeno_drivers, Adeno_samples)$test
  p_values_Adeno['SWEET_lung', j] <- test_hyper(SWEET_lung, Adeno_drivers, Adeno_samples)$test
  
  p_values_Squamous['LIONESS_lung', j] <- test_hyper(LIONESS_lung, Squamous_drivers, Squamous_samples)$test
  p_values_Squamous['SSN_lung', j] <- test_hyper(SSN_lung, Squamous_drivers, Squamous_samples)$test
  p_values_Squamous['ssPCC_lung', j] <- test_hyper(iENA_lung, Squamous_drivers, Squamous_samples)$test
  p_values_Squamous['CSN_lung', j] <- test_hyper(CSN_lung, Squamous_drivers, Squamous_samples)$test
  p_values_Squamous['SSPGI_lung', j] <- test_hyper(SSPGI_lung, Squamous_drivers, Squamous_samples)$test
  p_values_Squamous['SWEET_lung', j] <- test_hyper(SWEET_lung, Squamous_drivers, Squamous_samples)$test
  
  p_values_NSCLC['LIONESS_lung', j] <- test_hyper(LIONESS_lung, NSCLC_drivers, NSCLC_samples)$test
  p_values_NSCLC['SSN_lung', j] <- test_hyper(SSN_lung, NSCLC_drivers, NSCLC_samples)$test
  p_values_NSCLC['ssPCC_lung', j] <- test_hyper(iENA_lung, NSCLC_drivers, NSCLC_samples)$test
  p_values_NSCLC['CSN_lung', j] <- test_hyper(CSN_lung, NSCLC_drivers, NSCLC_samples)$test
  p_values_NSCLC['SSPGI_lung', j] <- test_hyper(SSPGI_lung, NSCLC_drivers, NSCLC_samples)$test
  p_values_NSCLC['SWEET_lung', j] <- test_hyper(SWEET_lung, NSCLC_drivers, NSCLC_samples)$test
  
  
  # Now test enrichment under a hypergeometric distribution for both SCLC and NSCLC
  fold_changes_SCLC['LIONESS_lung', j] <- test_hyper(LIONESS_lung, SCLC_drivers, SCLC_samples)$fold_change
  fold_changes_SCLC['SSN_lung', j] <- test_hyper(SSN_lung, SCLC_drivers, SCLC_samples)$fold_change
  fold_changes_SCLC['ssPCC_lung', j] <- test_hyper(iENA_lung, SCLC_drivers, SCLC_samples)$fold_change
  fold_changes_SCLC['CSN_lung', j] <- test_hyper(CSN_lung, SCLC_drivers, SCLC_samples)$fold_change
  fold_changes_SCLC['SSPGI_lung', j] <- test_hyper(SSPGI_lung, SCLC_drivers, SCLC_samples)$fold_change
  fold_changes_SCLC['SWEET_lung', j] <- test_hyper(SWEET_lung, SCLC_drivers, SCLC_samples)$fold_change
  
  fold_changes_Adeno['LIONESS_lung', j] <- test_hyper(LIONESS_lung, Adeno_drivers, Adeno_samples)$fold_change
  fold_changes_Adeno['SSN_lung', j] <- test_hyper(SSN_lung, Adeno_drivers, Adeno_samples)$fold_change
  fold_changes_Adeno['ssPCC_lung', j] <- test_hyper(iENA_lung, Adeno_drivers, Adeno_samples)$fold_change
  fold_changes_Adeno['CSN_lung', j] <- test_hyper(CSN_lung, Adeno_drivers, Adeno_samples)$fold_change
  fold_changes_Adeno['SSPGI_lung', j] <- test_hyper(SSPGI_lung, Adeno_drivers, Adeno_samples)$fold_change
  fold_changes_Adeno['SWEET_lung', j] <- test_hyper(SWEET_lung, Adeno_drivers, Adeno_samples)$fold_change
  
  fold_changes_Squamous['LIONESS_lung', j] <- test_hyper(LIONESS_lung, Squamous_drivers, Squamous_samples)$fold_change
  fold_changes_Squamous['SSN_lung', j] <- test_hyper(SSN_lung, Squamous_drivers, Squamous_samples)$fold_change
  fold_changes_Squamous['ssPCC_lung', j] <- test_hyper(iENA_lung, Squamous_drivers, Squamous_samples)$fold_change
  fold_changes_Squamous['CSN_lung', j] <- test_hyper(CSN_lung, Squamous_drivers, Squamous_samples)$fold_change
  fold_changes_Squamous['SSPGI_lung', j] <- test_hyper(SSPGI_lung, Squamous_drivers, Squamous_samples)$fold_change
  fold_changes_Squamous['SWEET_lung', j] <- test_hyper(SWEET_lung, Squamous_drivers, Squamous_samples)$fold_change
  
  fold_changes_NSCLC['LIONESS_lung', j] <- test_hyper(LIONESS_lung, NSCLC_drivers, NSCLC_samples)$fold_change
  fold_changes_NSCLC['SSN_lung', j] <- test_hyper(SSN_lung, NSCLC_drivers, NSCLC_samples)$fold_change
  fold_changes_NSCLC['ssPCC_lung', j] <- test_hyper(iENA_lung, NSCLC_drivers, NSCLC_samples)$fold_change
  fold_changes_NSCLC['CSN_lung', j] <- test_hyper(CSN_lung, NSCLC_drivers, NSCLC_samples)$fold_change
  fold_changes_NSCLC['SSPGI_lung', j] <- test_hyper(SSPGI_lung, NSCLC_drivers, NSCLC_samples)$fold_change
  fold_changes_NSCLC['SWEET_lung', j] <- test_hyper(SWEET_lung, NSCLC_drivers, NSCLC_samples)$fold_change
  
}

save.image('../results/Rebuttal/top200_connected_nodes/results_top200_lung_Intogen_drivers_withSWEET.RData')
load('../results/top200_connected_nodes/lung/results_top200.RData')

names(SSN_lung) <- gsub('X','', names(SSN_lung))
names(LIONESS_lung) <- gsub('X','', names(LIONESS_lung))
names(iENA_lung) <- gsub('X','', names(iENA_lung))
names(CSN_lung) <- gsub('X','', names(CSN_lung))
names(SSPGI_lung) <- gsub('X','', names(SSPGI_lung))
names(SWEET_lung) <- gsub('X','', names(SWEET_lung))

# Define your sample names
sample_names <- names(LIONESS_lung)

# Create an empty data frame to store the overlap results
overlap_results <- list()

# Loop through each sample and calculate the overlap
for (sample in sample_names) {
  
  #sample <- '0010'
  
  # Get the sets of unique values for each method
  lioness_set <- unlist(LIONESS_lung[[sample]])
  ssn_set <- unlist(SSN_lung[[sample]])
  sweet_set <- unlist(SWEET_lung[[sample]])
  iena_set <- unlist(iENA_lung[[sample]])
  csn_set <- unlist(CSN_lung[[sample]])
  sspgi_set <- unlist(SSPGI_lung[[sample]])
  
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

LIONESS_intersection <- Reduce(intersect, unname(LIONESS_lung)) # zero genes overlapping across all samples
SSN_intersection <- Reduce(intersect, SSN_lung) # zero genes
CSN_intersection <- Reduce(intersect, CSN_lung) # 96 genes (so really a lot if you think about the fact that we use top200 as hubs)
iENA_intersection <- Reduce(intersect, iENA_lung) # zero genes
SSPGI_intersection <- Reduce(intersect, SSPGI_lung) # 18 genes

# Is there more overlap between samples of a single subtype?
overlap_within_subtype <- function(top200list){
  
  names(top200list) <- gsub('X', '', names(top200list))
  #top200list <- LIONESS_lung
  top200_NSCLC <- top200list[names(top200list) %in% NSCLC_samples$DepMap_ID]
  top200_SCLC <- top200list[names(top200list) %in% SCLC_samples$DepMap_ID]
  #print(length(Reduce(intersect, unname(top200_NSCLC))))
  #print(length(Reduce(intersect, unname(top200_SCLC))))
  print('SCLC')
  print(length(unique(unlist(top200_SCLC))))
  print('NSCLC')
  print(length(unique(unlist(top200_NSCLC))))
  # print(paste0(length(unique(unlist(top200_NSCLC))), ' hub genes in NSCLC for ', top200list))
  # print(paste0(length(unique(unlist(top200_SCLC))), ' hub genes in SCLC for ', top200list))
}

overlap_within_subtype(SSN_lung)
overlap_within_subtype(LIONESS_lung)
overlap_within_subtype(iENA_lung)
overlap_within_subtype(CSN_lung)
overlap_within_subtype(SSPGI_lung)

# load('../results/top200_connected_nodes/lung/results5_500.RData')
#p_values_NSCLC <- fread('../results/top200_connected_nodes/lung/p_values_NSCLC.txt', header=FALSE)
#p_values_SCLC <- fread('../results/top200_connected_nodes/lung/p_values_SCLC.txt', header=FALSE)
# Create plot of p-values per subtype
rownames(p_values_NSCLC) <- c('LIONESS', 'SSN', 'ssPCC', 'CSN', 'SSPGI')
NSCLC_t <- as.data.frame(p_values_NSCLC)
NSCLC_t$method <-row.names(NSCLC_t)
NSCLC_t <- melt(NSCLC_t); colnames(NSCLC_t) <- c('Method', 'Number of hubs per sample', 'p-value')

pdf('../results/top200_connected_nodes/lung/NSCLC_p_values.pdf')
ggplot(NSCLC_t, aes(x=`Number of hubs per sample`, y=`p-value`, color = Method, group = Method)) + geom_line() + theme_bw() +
  scale_x_discrete(breaks=seq(1,200,1)) + geom_hline(yintercept=0.05) + scale_y_continuous(breaks=c(0.05, 0.25, 0.5, 0.75, 1))
dev.off()

rownames(p_values_SCLC) <- c('LIONESS', 'SSN', 'ssPCC', 'CSN', 'SSPGI')
SCLC_t <- as.data.frame(p_values_SCLC)
SCLC_t$method <-row.names(SCLC_t)
SCLC_t <- melt(SCLC_t); colnames(SCLC_t) <- c('Method', 'Number of hubs per sample', 'p-value')

pdf('../results/top200_connected_nodes/lung/SCLC_p_values.pdf')
ggplot(SCLC_t, aes(x=`Number of hubs per sample`, y=`p-value`, color = Method, group = Method)) + geom_line() + theme_bw() +
  scale_x_discrete(breaks=seq(5,500,50)) + geom_hline(yintercept=0.05) + scale_y_continuous(breaks=c(0.05, 0.25, 0.5, 0.75, 1))
dev.off()


          





