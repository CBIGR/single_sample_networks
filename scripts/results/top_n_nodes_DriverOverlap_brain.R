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
p_values_HGG <- data.frame(matrix(nrow = 5, ncol=length(number_of_hubs)))
rownames(p_values_HGG) <- c('LIONESS_brain', 'SSN_brain', 'ssPCC_brain', 'CSN_brain', 'SSPGI_brain'); colnames(p_values_HGG) <- number_of_hubs
p_values_MBL <- data.frame(matrix(nrow = 5, ncol=length(number_of_hubs)))
rownames(p_values_MBL) <- c('LIONESS_brain', 'SSN_brain', 'ssPCC_brain', 'CSN_brain', 'SSPGI_brain'); colnames(p_values_MBL) <- number_of_hubs

fold_changes_HGG <- data.frame(matrix(nrow = 5, ncol=length(number_of_hubs)))
rownames(fold_changes_HGG) <- c('LIONESS_brain', 'SSN_brain', 'ssPCC_brain', 'CSN_brain', 'SSPGI_brain'); colnames(fold_changes_HGG) <- number_of_hubs
fold_changes_MBL <- data.frame(matrix(nrow = 5, ncol=length(number_of_hubs)))
rownames(fold_changes_MBL) <- c('LIONESS_brain', 'SSN_brain', 'ssPCC_brain', 'CSN_brain', 'SSPGI_brain'); colnames(fold_changes_MBL) <- number_of_hubs


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
  # save.image(paste0('../results/top200_connected_nodes/top',k,'nodes_brain.RData'))
  
  overlap_df <- data.frame(matrix(NA, nrow = length(names(LIONESS_brain)), ncol = 28), row.names = names(LIONESS_brain))
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
  
  for (i in 1:67){ # 67 samples in brain
    # i <- 5
    # LIONESS <- convertToSymbol(LIONESS_brain[[i]])
    # iENA <- convertToSymbol(iENA_brain[[i]])
    # SSPGI <- convertToSymbol(SSPGI_brain[[i]])
    # SSN <- convertToSymbol(SSN_brain[[i]])
    # CSN <- convertToSymbol(CSN_brain[[i]])
    
    listInput <- list(LIONESS = LIONESS_brain[[i]], iENA = iENA_brain[[i]], SSN = SSN_brain[[i]], SSPGI=SSPGI_brain[[i]], CSN=CSN_brain[[i]])
    plot <- upset(fromList(listInput), order.by = "freq")
    # print(plot)
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
  # pdf(paste0('../results/top200_connected_nodes/brain/Boxplots/Top_',k,'_connected_nodes_brain_overlap.pdf'))
  # print(ggplot(data_for_boxplot, aes(x = `Method(s)`, y=`Number of overlapping hubs`, fill = Group)) + geom_boxplot() +
  #         ggtitle('Overlap between the top 200 most connected nodes') + scale_colour_manual(values=cbbPalette) +
  #         theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=45, hjust=1, size = 6), 
  #               panel.background = element_rect(fill = 'white', colour = 'white'),
  #               # panel.grid.minor = element_line(size=0.2, linetype='solid', colour='black'), 
  #               panel.grid.major = element_line(size=0.1, linetype='solid', colour='grey'),
  #               axis.line = element_line(size = 0.5, linetype = 'solid', colour = 'black'))) 
  # dev.off()
  # 
  # Read in sample metadata and group samples based on subtype
  sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,19,24)]
  sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  MBL_samples <- sample_metadata[sample_metadata$Subtype == 'Medulloblastoma', ]
  HGG_samples <- sample_metadata[sample_metadata$Subtype == 'Glioblastoma', ]
  
   #Providing cancer type as a variable doesn't work with paste0 so change this manually
  type <- 'Medulloblastoma'; print(type)
  LIONESS_MBL <- overlapWithDrivers(LIONESS_brain, type, MBL_drivers)
  SSN_MBL <- overlapWithDrivers(SSN_brain, type, MBL_drivers)
  CSN_MBL <- overlapWithDrivers(CSN_brain, type, MBL_drivers)
  ssPCC_MBL <- overlapWithDrivers(iENA_brain, type, MBL_drivers)
  SSPGI_MBL <- overlapWithDrivers(SSPGI_brain, type, MBL_drivers)
  # pdf(paste0('../results/top200_connected_nodes/brain/Venn/',type,'_top', k, 'nodes_driverOverlap.pdf'))
  # print(ggarrange(SSN_MBL, LIONESS_MBL, ssPCC_MBL, LIONESS_MBL, SSPGI_MBL, nrow=3, ncol=2, labels = c('SSN', 'LIONESS', 'iENA', 'CSN', 'SSPGI'), font.label = list(size = 8)))
  # dev.off()
  
  type <- 'Glioblastoma'; print(type)
  LIONESS_HGG <- overlapWithDrivers(LIONESS_brain, type, HGG_drivers)
  SSN_HGG <- overlapWithDrivers(SSN_brain, type, HGG_drivers)
  CSN_HGG <- overlapWithDrivers(CSN_brain, type, HGG_drivers)
  ssPCC_HGG <- overlapWithDrivers(iENA_brain, type, HGG_drivers)
  SSPGI_HGG <- overlapWithDrivers(SSPGI_brain, type, HGG_drivers)
  # pdf(paste0('../results/top200_connected_nodes/brain/Venn/',type,'_top', k, 'nodes_driverOverlap.pdf'))
  # print(ggarrange(SSN_HGG, LIONESS_HGG, ssPCC_HGG, CSN_HGG, SSPGI_HGG, nrow=3, ncol=2, labels = c('SSN', 'LIONESS', 'iENA', 'CSN', 'SSPGI'), font.label = list(size = 8)))
  # dev.off()
  
  p_values_MBL['LIONESS_brain', j] <- test_hyper(LIONESS_brain, MBL_drivers, MBL_samples)$test
  p_values_MBL['SSN_brain', j] <- test_hyper(SSN_brain,  MBL_drivers, MBL_samples)$test
  p_values_MBL['ssPCC_brain', j] <- test_hyper(iENA_brain, MBL_drivers, MBL_samples)$test
  p_values_MBL['CSN_brain', j] <- test_hyper(CSN_brain, MBL_drivers, MBL_samples)$test
  p_values_MBL['SSPGI_brain', j] <- test_hyper(SSPGI_brain, MBL_drivers, MBL_samples)$test
  
  p_values_HGG['LIONESS_brain', j] <- test_hyper(LIONESS_brain, HGG_drivers, HGG_samples)$test
  p_values_HGG['SSN_brain', j] <- test_hyper(SSN_brain, HGG_drivers, HGG_samples)$test
  p_values_HGG['ssPCC_brain', j] <- test_hyper(iENA_brain, HGG_drivers, HGG_samples)$test
  p_values_HGG['CSN_brain', j] <- test_hyper(CSN_brain, HGG_drivers, HGG_samples)$test
  p_values_HGG['SSPGI_brain', j] <- test_hyper(SSPGI_brain, HGG_drivers, HGG_samples)$test
  
  
  fold_changes_MBL['LIONESS_brain', j] <- test_hyper(LIONESS_brain, MBL_drivers, MBL_samples)$fold_change
  fold_changes_MBL['SSN_brain', j] <- test_hyper(SSN_brain,  MBL_drivers, MBL_samples)$fold_change
  fold_changes_MBL['ssPCC_brain', j] <- test_hyper(iENA_brain, MBL_drivers, MBL_samples)$fold_change
  fold_changes_MBL['CSN_brain', j] <- test_hyper(CSN_brain, MBL_drivers, MBL_samples)$fold_change
  fold_changes_MBL['SSPGI_brain', j] <- test_hyper(SSPGI_brain, MBL_drivers, MBL_samples)$fold_change
  
  fold_changes_HGG['LIONESS_brain', j] <- test_hyper(LIONESS_brain, HGG_drivers, HGG_samples)$fold_change
  fold_changes_HGG['SSN_brain', j] <- test_hyper(SSN_brain, HGG_drivers, HGG_samples)$fold_change
  fold_changes_HGG['ssPCC_brain', j] <- test_hyper(iENA_brain, HGG_drivers, HGG_samples)$fold_change
  fold_changes_HGG['CSN_brain', j] <- test_hyper(CSN_brain, HGG_drivers, HGG_samples)$fold_change
  fold_changes_HGG['SSPGI_brain', j] <- test_hyper(SSPGI_brain, HGG_drivers, HGG_samples)$fold_change
  
  write.table(p_values_HGG, '../results/top200_connected_nodes/brain/p_values_HGG.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  write.table(p_values_MBL, '../results/top200_connected_nodes/brain/p_values_MBL.txt', row.names = TRUE, sep = '\t', quote = FALSE)

  write.table(fold_changes_HGG, '../results/top200_connected_nodes/brain/fold_changes_HGG.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  write.table(fold_changes_MBL, '../results/top200_connected_nodes/brain/fold_changes_MBL.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  
  
  # write.table(p_values_HGG, '../results/top200_connected_nodes/brain/p_values_HGG_vs_allDrivers.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  # write.table(p_values_MBL, '../results/top200_connected_nodes/brain/p_values_MBL_vs_allDrivers.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  # 
  # write.table(fold_changes_HGG, '../results/top200_connected_nodes/brain/fold_changes_HGG_vs_allDrivers.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  # write.table(fold_changes_MBL, '../results/top200_connected_nodes/brain/fold_changes_MBL_vs_allDrivers.txt', row.names = TRUE, sep = '\t', quote = FALSE)
  
}

save.image('../results/top200_connected_nodes/brain/results200.RData')
#load('../results/top200_connected_nodes/brain/results5_500.RData')

LIONESS_intersection <- Reduce(intersect, unname(LIONESS_brain)) # zero genes overlapping across all samples
SSN_intersection <- Reduce(intersect, SSN_brain) # zero genes
CSN_intersection <- Reduce(intersect, CSN_brain) # 76 genes (so really a lot if you think about the fact that we use top200 as hubs)
iENA_intersection <- Reduce(intersect, iENA_brain) # zero genes
SSPGI_intersection <- Reduce(intersect, SSPGI_brain) # 21 genes

overlap_within_subtype <- function(top200list){
  # top200list <- LIONESS_brain
  top200_MBL <- top200list[names(top200list) %in% MBL_samples$DepMap_ID]
  top200_HGG <- top200list[names(top200list) %in% HGG_samples$DepMap_ID]
  print(length(Reduce(intersect, unname(top200_MBL))))
  print(length(Reduce(intersect, unname(top200_HGG))))
}

overlap_within_subtype(SSN_brain)
overlap_within_subtype(LIONESS_brain)
overlap_within_subtype(iENA_brain)
overlap_within_subtype(CSN_brain)
overlap_within_subtype(SSPGI_brain)

# # Create plot of p-values per subtype
rownames(p_values_HGG) <- c('LIONESS', 'SSN', 'ssPCC', 'CSN', 'SSPGI')
HGG_t <- as.data.frame(p_values_HGG)
HGG_t$method <-row.names(HGG_t)
HGG_t <- melt(HGG_t); colnames(HGG_t) <- c('Method', 'Number of hubs per sample', 'p-value')

pdf('../results/top200_connected_nodes/brain/HGG_p_values.pdf')
ggplot(HGG_t, aes(x=`Number of hubs per sample`, y=`p-value`, color = Method, group = Method)) + geom_line() + theme_bw() +
  scale_x_discrete(breaks=seq(0,200,10)) + geom_hline(yintercept=0.05) + scale_y_continuous(breaks=c(0.05, 0.25, 0.5, 0.75,1))
dev.off()

rownames(p_values_MBL) <- c('LIONESS', 'SSN', 'ssPCC', 'CSN', 'SSPGI')
MBL_t <- as.data.frame(p_values_MBL)
MBL_t$method <-row.names(MBL_t)
MBL_t <- melt(MBL_t); colnames(MBL_t) <- c('Method', 'Number of hubs per sample', 'p-value')

pdf('../results/top200_connected_nodes/brain/MBL_p_values.pdf')
ggplot(MBL_t, aes(x=`Number of hubs per sample`, y=`p-value`, color = Method, group = Method)) + geom_line() + theme_bw() +
  scale_x_discrete(breaks=seq(0,200,10)) + geom_hline(yintercept=0.05) + scale_y_continuous(breaks=c(0.05, 0.25, 0.5, 0.75,1))
dev.off()

