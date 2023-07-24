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
  # network <- paste0(base, "SSPGI/SSPGI_lung_edges_2_75_HumanNet_top_25000_edges.csv")
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
  
  # top200list <- LIONESS_lung
  # driverlist <- NSCLC_drivers
  # cancerType <- 'NSCLC'
  
  names(top200list) <- gsub('X', '', names(top200list))
  for (i in 1:length(names(top200list))){ # Convert top200 nodes to gene symbols
    top200list[[i]] <- convertToSymbol(top200list[[i]])
  }
  
  # Select samples belonging to the subtupe of interest, extract union of top200 nodes
  # sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,24)]
  # sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  samples <- sample_metadata[sample_metadata$lineage_subtype == cancerType, ]
  top200list <- top200list[names(top200list) %in% samples$DepMap_ID]
  nodes <- unique(unlist(unname(top200list))) #Union of all driver genes found in samples of the subtype of interest
  
  
  if (cancerType == 'SCLC' | cancerType == 'NSCLC'){
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
  # method <- SSN_lung
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
NSCLC_drivers_intogen <- fread('./IntOGen-DriverGenes_NSCLC.tsv')[,1]
SCLC_drivers_census <- fread('./Census_SCLC_19_06_2023.tsv')[,1]
NSCLC_drivers_census <- fread('./Census_NSCLC_19_06_2023.tsv')[,1]

SCLC_drivers <- rbind(SCLC_drivers_intogen, SCLC_drivers_census, use.names=FALSE)
NSCLC_drivers <- rbind(NSCLC_drivers_intogen, NSCLC_drivers_census, use.names=FALSE)

all_drivers_intogen <- fread('./IntOGen-DriverGenes_all.tsv')[,1]
all_drivers_census <- fread('./Census_all_23_06_2023.tsv')[,1]
all_drivers <- rbind(all_drivers_intogen, all_drivers_census, use.names=FALSE)


HumanNet_genes <- unique(as.vector(as.matrix(aggregate[,c(1,2)])))
tmp1 <- intersect(HumanNet_genes, NSCLC_drivers$Symbol)
tmp2 <- intersect(HumanNet_genes, SCLC_drivers$Symbol)
length(union(tmp1, tmp2))


# number_of_hubs <- seq(5,500,10)
#number_of_hubs <- seq(1,200,5)
number_of_hubs <- c(200)
p_values_SCLC <- data.frame(matrix(nrow = 5, ncol=length(number_of_hubs)))
rownames(p_values_SCLC) <- c('LIONESS_lung', 'SSN_lung', 'ssPCC_lung', 'CSN_lung', 'SSPGI_lung'); colnames(p_values_SCLC) <- number_of_hubs
p_values_NSCLC <- data.frame(matrix(nrow = 5, ncol=length(number_of_hubs)))
rownames(p_values_NSCLC) <- c('LIONESS_lung', 'SSN_lung', 'ssPCC_lung', 'CSN_lung', 'SSPGI_lung'); colnames(p_values_NSCLC) <- number_of_hubs

fold_changes_SCLC <- data.frame(matrix(nrow = 5, ncol=length(number_of_hubs)))
rownames(fold_changes_SCLC) <- c('LIONESS_lung', 'SSN_lung', 'ssPCC_lung', 'CSN_lung', 'SSPGI_lung'); colnames(fold_changes_SCLC) <- number_of_hubs
fold_changes_NSCLC <- data.frame(matrix(nrow = 5, ncol=length(number_of_hubs)))
rownames(fold_changes_NSCLC) <- c('LIONESS_lung', 'SSN_lung', 'ssPCC_lung', 'CSN_lung', 'SSPGI_lung'); colnames(fold_changes_NSCLC) <- number_of_hubs

aggregate <- fread('./aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv', fill = TRUE, header = TRUE, data.table = FALSE) #Should we use the complete aggregate network or the top 25k edges?
aggregate$reg <- convertToSymbol(aggregate$reg)
aggregate$tar <- convertToSymbol(aggregate$tar)
HumanNet_genes <- unique(as.vector(as.matrix(aggregate[,c(1,2)])))


# number_of_hubs <- c(5, 200)
for (j in (1:length(number_of_hubs))){
  # j <- 1
  k <- number_of_hubs[j]
  # Find hub genes
  LIONESS_lung <- top200nodes(paste0(base, "LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv"),k, offset=TRUE)
  CSN_lung <- top200nodes(paste0(base, "CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv"), k, offset=FALSE)
  SSN_lung <- top200nodes(paste0(base,"SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), k, offset=TRUE)
  iENA_lung <- top200nodes(paste0(base,"ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), k, offset=TRUE)
  SSPGI_lung <- top200nodes(paste0(base,"SSPGI/SSPGI_lung_edges_2_75_HumanNet_top_25000_edges_symbol.csv"), k, offset=TRUE)
  # save.image(paste0('../results/top200_connected_nodes/lung/top',k,'nodes.RData'))
  
  overlap_df <- data.frame(matrix(NA, nrow = length(names(LIONESS_lung)), ncol = 28), row.names = names(LIONESS_lung))
  colnames(overlap_df) <- c('LIONESS_SSN_CSN_SSPGI_iENA',
                            'LIONESS_SSN_CSN_SSPGI', # Options in case overlap size is 4 groups
                            'LIONESS_SSN_CSN_iENA',
                            'LIONESS_SSN_SSPGI_iENA',
                            'LIONESS_CSN_SSPGI_iENA',
                            'SSN_CSN_SSPGI_iENA',
                            'LIONESS_SSN_CSN', # Options in case overlap size is 3
                            'LIONESS_CSN_SSPGI', 
                            'LIONESS_SSPGI_iENA',
                            'LIONESS_CSN_iENA', # to add
                            'SSN_CSN_SSPGI',
                            'SSN_SSPGI_iENA',
                            'CSN_SSPGI_iENA', 
                            'CSN_SSN_iENA', # Start here
                            'LIONESS_CSN', 'LIONESS_SSN', 'LIONESS_SSPGI', 'LIONESS_iENA', # 2 matches
                            'SSN_CSN', 'CSN_SSPGI', 'CSN_iENA',
                            'SSN_SSPGI', 'SSN_iENA',
                            'LIONESS', 'CSN', 'SSN', 'SSPGI', 'iENA' # Only one match
  )
  
  
  for (i in 1:86){ # 86 samples in lung
    # LIONESS <- convertToSymbol(LIONESS_lung[[i]])
    # iENA <- convertToSymbol(iENA_lung[[i]])
    # SSPGI <- convertToSymbol(SSPGI_lung[[i]])
    # SSN <- convertToSymbol(SSN_lung[[i]])
    # CSN <- convertToSymbol(CSN_lung[[i]])
    
    listInput <- list(LIONESS = LIONESS_lung[[i]], iENA = iENA_lung[[i]], SSN = SSN_lung[[i]], SSPGI=SSPGI_lung[[i]], CSN=CSN_lung[[i]])
    plot <- upset(fromList(listInput), order.by = "freq")
    plot
    # From this upSet plot you can extract a binary list with overlaps
    overlap_table <- plot$New_data
    overlap_table <- overlap_table[,c(1,3,5,4,2)]
    # I can add rownames to these tables in case I would need them later on, but skip for now
    # Now calculate the number of overlapping features for each column in the overview dataframe
    overlap_df[i, 'LIONESS_SSN_CSN_SSPGI_iENA'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 1 & overlap_table[,3] == 1 & overlap_table[,4] == 1 & overlap_table[,5] == 1, ])
    overlap_df[i, 'LIONESS_SSN_CSN_SSPGI'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 1 & overlap_table[,3] == 1 & overlap_table[,4] == 1 & overlap_table[,5] == 0, ])
    overlap_df[i, 'LIONESS_SSN_CSN_iENA'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 1 & overlap_table[,3] == 1 & overlap_table[,4] == 0 & overlap_table[,5] == 1, ])
    overlap_df[i, 'LIONESS_SSN_SSPGI_iENA'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 1 & overlap_table[,3] == 0 & overlap_table[,4] == 1 & overlap_table[,5] == 1, ])
    overlap_df[i, 'LIONESS_CSN_SSPGI_iENA'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 0 & overlap_table[,3] == 1 & overlap_table[,4] == 1 & overlap_table[,5] == 1, ])
    overlap_df[i, 'SSN_CSN_SSPGI_iENA'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 1 & overlap_table[,3] == 1 & overlap_table[,4] == 1 & overlap_table[,5] == 1, ])
    
    overlap_df[i, 'LIONESS_SSN_CSN'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 1 & overlap_table[,3] == 1 & overlap_table[,4] == 0 & overlap_table[,5] == 0, ])
    overlap_df[i, 'LIONESS_CSN_SSPGI'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 0 & overlap_table[,3] == 1 & overlap_table[,4] == 1 & overlap_table[,5] == 0, ])
    overlap_df[i, 'LIONESS_SSN_iENA'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 1 & overlap_table[,3] == 0 & overlap_table[,4] == 0 & overlap_table[,5] == 1, ])
    overlap_df[i, 'LIONESS_CSN_iENA'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 0 & overlap_table[,3] == 1 & overlap_table[,4] == 0 & overlap_table[,5] == 0, ])
    overlap_df[i, 'LIONESS_SSPGI_iENA'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 0 & overlap_table[,3] == 0 & overlap_table[,4] == 1 & overlap_table[,5] == 1, ])
    overlap_df[i, 'LIONESS_CSN_iENA'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 1 & overlap_table[,3] == 0 & overlap_table[,4] == 0 & overlap_table[,5] == 1, ])
    overlap_df[i, 'SSN_CSN_SSPGI'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 1 & overlap_table[,3] == 1 & overlap_table[,4] == 1 & overlap_table[,5] == 0, ])
    overlap_df[i, 'SSN_SSPGI_iENA'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 1 & overlap_table[,3] == 0 & overlap_table[,4] == 1 & overlap_table[,5] == 1, ])
    overlap_df[i, 'CSN_SSPGI_iENA'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 0 & overlap_table[,3] == 1 & overlap_table[,4] == 1 & overlap_table[,5] == 1, ])
    overlap_df[i, 'CSN_SSN_iENA'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 1 & overlap_table[,3] == 1 & overlap_table[,4] == 0 & overlap_table[,5] == 1, ])
    
    overlap_df[i, 'LIONESS_SSN'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 1 & overlap_table[,3] == 0 & overlap_table[,4] == 0 & overlap_table[,5] == 0, ])
    overlap_df[i, 'LIONESS_CSN'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 0 & overlap_table[,3] == 1 & overlap_table[,4] == 0 & overlap_table[,5] == 0, ])
    overlap_df[i, 'LIONESS_SSPGI'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 0 & overlap_table[,3] == 0 & overlap_table[,4] == 1 & overlap_table[,5] == 0, ])
    overlap_df[i, 'LIONESS_iENA'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 0 & overlap_table[,3] == 0 & overlap_table[,4] == 0 & overlap_table[,5] == 1, ])
    overlap_df[i, 'SSN_CSN'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 1 & overlap_table[,3] == 1 & overlap_table[,4] == 0 & overlap_table[,5] == 0, ])
    overlap_df[i, 'SSN_SSPGI'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 1 & overlap_table[,3] == 0 & overlap_table[,4] == 1 & overlap_table[,5] == 0, ])
    overlap_df[i, 'SSN_iENA'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 1 & overlap_table[,3] == 0 & overlap_table[,4] == 0 & overlap_table[,5] == 1, ])
    overlap_df[i, 'CSN_SSPGI'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 0 & overlap_table[,3] == 1 & overlap_table[,4] == 1 & overlap_table[,5] == 0, ])
    overlap_df[i, 'CSN_iENA'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 0 & overlap_table[,3] == 1 & overlap_table[,4] == 0 & overlap_table[,5] == 1, ])
    overlap_df[i, 'SSPGI_iENA'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 0 & overlap_table[,3] == 0 & overlap_table[,4] == 1 & overlap_table[,5] == 1, ])
    
    overlap_df[i, 'LIONESS'] <- nrow(overlap_table[overlap_table[,1] == 1 & overlap_table[,2] == 0 & overlap_table[,3] == 0 & overlap_table[,4] == 0 & overlap_table[,5] == 0, ])
    overlap_df[i, 'SSN'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 1 & overlap_table[,3] == 0 & overlap_table[,4] == 0 & overlap_table[,5] == 0, ])
    overlap_df[i, 'CSN'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 0 & overlap_table[,3] == 1 & overlap_table[,4] == 0 & overlap_table[,5] == 0, ])
    overlap_df[i, 'SSPGI'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 0 & overlap_table[,3] == 0 & overlap_table[,4] == 1 & overlap_table[,5] == 0, ])
    overlap_df[i, 'iENA'] <- nrow(overlap_table[overlap_table[,1] == 0 & overlap_table[,2] == 0 & overlap_table[,3] == 0 & overlap_table[,4] == 0 & overlap_table[,5] == 1, ])
    
  } # Create overlap dataframe
  
  filtered_overlap_df <- overlap_df[, colSums(overlap_df) > -1]
  filtered_overlap_df <- filtered_overlap_df[,c(24:28, 15:23, 30 ,29, 7:14, 2:6,1)]
  
  data_for_boxplot <- reshape2::melt(filtered_overlap_df)
  colnames(data_for_boxplot) <- c('Method(s)', 'Number of overlapping genes')
  data_for_boxplot$group <- NA
  for (i in (1:nrow(data_for_boxplot))){ # group boxplots according to number of overlaps
    data_for_boxplot[i, 'group'] <- length(str_split(data_for_boxplot[i, 'Method(s)'], '_')[[1]])
  }
  data_for_boxplot$group <- as.factor(data_for_boxplot$group)
  colnames(data_for_boxplot) <- c('Method(s)', 'Number of overlapping hubs', 'Group')
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # pdf(paste0('../results/top200_connected_nodes/lung/Boxplots/Top_',k,'_connected_nodes_lung_overlap.pdf'))
  # print(ggplot(data_for_boxplot, aes(x = `Method(s)`, y=`Number of overlapping hubs`, fill = Group)) + geom_boxplot() +
  #         ggtitle('Overlap between the top 200 most connected nodes') + scale_colour_manual(values=cbbPalette) +
  #         theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1, size = 6), 
  #               panel.background = element_rect(fill = 'white', colour = 'white'),
  #               # panel.grid.minor = element_line(size=0.2, linetype='solid', colour='black'), 
  #               panel.grid.major = element_line(size=0.1, linetype='solid', colour='grey'),
  #               axis.line = element_line(size = 0.5, linetype = 'solid', colour = 'black'))) 
  # dev.off()
  
  # Read in sample metadata and group samples based on subtype
  sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,24)]
  sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  NSCLC_samples <- sample_metadata[sample_metadata$lineage_subtype == 'NSCLC', ]
  SCLC_samples <- sample_metadata[sample_metadata$lineage_subtype == 'SCLC', ]
  
  # Once for NSCLC and once for SCLC
  type <- 'SCLC'
  print('SCLC')
  LIONESS_SCLC <- overlapWithDrivers(LIONESS_lung, type, SCLC_drivers)
  SSN_SCLC <- overlapWithDrivers(SSN_lung, type, SCLC_drivers)
  CSN_SCLC <- overlapWithDrivers(CSN_lung, type, SCLC_drivers)
  ssPCC_SCLC <- overlapWithDrivers(iENA_lung, type, SCLC_drivers)
  SSPGI_SCLC <- overlapWithDrivers(SSPGI_lung, type, SCLC_drivers)
  
  # pdf(paste0('../results/top200_connected_nodes/lung/Venn/',type,'_top', k, 'nodes_driverOverlap.pdf'))
  # print(ggarrange(SSN_SCLC, LIONESS_SCLC, ssPCC_SCLC, CSN_SCLC, SSPGI_SCLC, nrow=3, ncol=2, labels = c('SSN', 'LIONESS', 'iENA', 'CSN', 'SSPGI'), font.label = list(size = 8)))
  # dev.off()
  
  # Once for NSCLC and once for SCLC
  type <- 'NSCLC'
  print('NSCLC')
  LIONESS_NSCLC <- overlapWithDrivers(LIONESS_lung, type, NSCLC_drivers)
  SSN_NSCLC <- overlapWithDrivers(SSN_lung, type, NSCLC_drivers)
  CSN_NSCLC <- overlapWithDrivers(CSN_lung, type, NSCLC_drivers)
  ssPCC_NSCLC <- overlapWithDrivers(iENA_lung, type, NSCLC_drivers)
  SSPGI_NSCLC <- overlapWithDrivers(SSPGI_lung, type, NSCLC_drivers)
  
  # pdf(paste0('../results/top200_connected_nodes/lung/Venn/',type,'_top', k, 'nodes_driverOverlap.pdf'))
  # print(ggarrange(SSN_NSCLC, LIONESS_NSCLC, ssPCC_NSCLC, CSN_NSCLC, SSPGI_NSCLC, nrow=3, ncol=2, labels = c('SSN', 'LIONESS', 'iENA', 'CSN', 'SSPGI'), font.label = list(size = 8)))
  # dev.off()
  
  # Now test enrichment under a hypergeometric distribution for both SCLC and NSCLC
  p_values_SCLC['LIONESS_lung', j] <- test_hyper(LIONESS_lung, SCLC_drivers, SCLC_samples)$test
  p_values_SCLC['SSN_lung', j] <- test_hyper(SSN_lung, SCLC_drivers, SCLC_samples)$test
  p_values_SCLC['ssPCC_lung', j] <- test_hyper(iENA_lung, SCLC_drivers, SCLC_samples)$test
  p_values_SCLC['CSN_lung', j] <- test_hyper(CSN_lung, SCLC_drivers, SCLC_samples)$test
  p_values_SCLC['SSPGI_lung', j] <- test_hyper(SSPGI_lung, SCLC_drivers, SCLC_samples)$test
  
  p_values_NSCLC['LIONESS_lung', j] <- test_hyper(LIONESS_lung, NSCLC_drivers, NSCLC_samples)$test
  p_values_NSCLC['SSN_lung', j] <- test_hyper(SSN_lung, NSCLC_drivers, NSCLC_samples)$test
  p_values_NSCLC['ssPCC_lung', j] <- test_hyper(iENA_lung, NSCLC_drivers, NSCLC_samples)$test
  p_values_NSCLC['CSN_lung', j] <- test_hyper(CSN_lung, NSCLC_drivers, NSCLC_samples)$test
  p_values_NSCLC['SSPGI_lung', j] <- test_hyper(SSPGI_lung, NSCLC_drivers, NSCLC_samples)$test
  
  
  # Now test enrichment under a hypergeometric distribution for both SCLC and NSCLC and report fold change, just do the test twice, doesn't take long
  fold_changes_SCLC['LIONESS_lung', j] <- test_hyper(LIONESS_lung, SCLC_drivers, SCLC_samples)$fold_change
  fold_changes_SCLC['SSN_lung', j] <- test_hyper(SSN_lung, SCLC_drivers, SCLC_samples)$fold_change
  fold_changes_SCLC['ssPCC_lung', j] <- test_hyper(iENA_lung, SCLC_drivers, SCLC_samples)$fold_change
  fold_changes_SCLC['CSN_lung', j] <- test_hyper(CSN_lung, SCLC_drivers, SCLC_samples)$fold_change
  fold_changes_SCLC['SSPGI_lung', j] <- test_hyper(SSPGI_lung, SCLC_drivers, SCLC_samples)$fold_change
  
  fold_changes_NSCLC['LIONESS_lung', j] <- test_hyper(LIONESS_lung, NSCLC_drivers, NSCLC_samples)$fold_change
  fold_changes_NSCLC['SSN_lung', j] <- test_hyper(SSN_lung, NSCLC_drivers, NSCLC_samples)$fold_change
  fold_changes_NSCLC['ssPCC_lung', j] <- test_hyper(iENA_lung, NSCLC_drivers, NSCLC_samples)$fold_change
  fold_changes_NSCLC['CSN_lung', j] <- test_hyper(CSN_lung, NSCLC_drivers, NSCLC_samples)$fold_change
  fold_changes_NSCLC['SSPGI_lung', j] <- test_hyper(SSPGI_lung, NSCLC_drivers, NSCLC_samples)$fold_change
  
  write.table(p_values_SCLC, '../results/top200_connected_nodes/lung/p_values_SCLC.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  write.table(p_values_NSCLC, '../results/top200_connected_nodes/lung/p_values_NSCLC.txt', row.names = TRUE, sep = '\t', quote = FALSE)

  write.table(fold_changes_SCLC, '../results/top200_connected_nodes/lung/fold_changes_SCLC.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  write.table(fold_changes_NSCLC, '../results/top200_connected_nodes/lung/fold_changes_NSCLC.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  
  # write.table(p_values_SCLC, '../results/top200_connected_nodes/lung/p_values_SCLC_vs_allDrivers.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  # write.table(p_values_NSCLC, '../results/top200_connected_nodes/lung/p_values_NSCLC_vs_allDrivers.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  # 
  # write.table(fold_changes_SCLC, '../results/top200_connected_nodes/lung/fold_changes_SCLC_vs_allDrivers.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  # write.table(fold_changes_NSCLC, '../results/top200_connected_nodes/lung/fold_changes_NSCLC_vs_allDrivers.txt', row.names = TRUE, sep = '\t', quote = FALSE)
}

save.image('../results/top200_connected_nodes/lung/results_top200.RData')
load('../results/top200_connected_nodes/lung/results_top200.RData')

# we find 35 overlapping genes in the top200 nodes, which are these?

LIONESS_intersection <- Reduce(intersect, unname(LIONESS_lung)) # zero genes overlapping across all samples
SSN_intersection <- Reduce(intersect, SSN_lung) # zero genes
CSN_intersection <- Reduce(intersect, CSN_lung) # 96 genes (so really a lot if you think about the fact that we use top200 as hubs)
iENA_intersection <- Reduce(intersect, iENA_lung) # zero genes
SSPGI_intersection <- Reduce(intersect, SSPGI_lung) # 18 genes

# Is there more overlap between samples of a single subtype?
overlap_within_subtype <- function(top200list){
  # top200list <- LIONESS_lung
  top200_NSCLC <- top200list[names(top200list) %in% NSCLC_samples$DepMap_ID]
  top200_SCLC <- top200list[names(top200list) %in% SCLC_samples$DepMap_ID]
  print(length(Reduce(intersect, unname(top200_NSCLC))))
  print(length(Reduce(intersect, unname(top200_SCLC))))
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


          





