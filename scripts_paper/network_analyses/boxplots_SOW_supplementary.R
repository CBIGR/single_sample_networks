#!/usr/bin/Rscript

library(data.table)
library(dplyr)
library(limma)
library(stringr)
library(ggplot2)
library(tidyr)
library(EnhancedVolcano)
library(biomaRt)
library(edgeR)
library(ggpubr)
library(parallel)
library(tidyverse)
library(reshape2)
library(purrr)


setwd('/home/boris/Documents/PhD/single_sample/networks')

calculate_SOW <- function(networkfile, offset=FALSE){
  network <-  network <- fread(networkfile, header=TRUE, data.table=FALSE, fill=TRUE)
  colnames(network) <- gsub('X', '', colnames(network))
  if (offset == TRUE){
    network[,1] <- NULL
  }
  
  # Create dataframe to store sum of weights
  genes <- unique(as.vector(as.matrix(network[,c(1,2)])))
  SOW <- data.frame(matrix(nrow = length(genes), ncol = ncol(network)-2))
  rownames(SOW) <- genes
  colnames(SOW) <- colnames(network[, -c(1,2)])
  
  # Calculate sum of absolute edge weights for all genes in every sample
  for (i in 1:length(genes)){
    # i <- 1
    gene <- rownames(SOW)[i]
    ind <- which(network[,1] %in% gene)
    ind <- c(ind, which(network[,2] %in% gene))
    SOW[i,] <- colSums(abs(network[ind, -c(1,2)])) # This is the approach used in the basic figure now included in the paper, comment lines 36:42, 55 out to change
    # SOW[i,] <- colSums(network[ind, -c(1,2)])
  }
  rownames(SOW) <- str_replace(row.names(SOW), "\\s.*", "")
  rownames(SOW) <- str_replace(row.names(SOW), "\\(.*", "")
  return(SOW)
}

LIONESS_SOW <- calculate_SOW('./LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset=FALSE)
SSN_SOW <- calculate_SOW('./SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset=FALSE)
CSN_SOW <-calculate_SOW('./CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset=FALSE)
iENA_SOW <-calculate_SOW('./ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv',offset=FALSE)
SSPGI_SOW <-calculate_SOW('./SSPGI/SSPGI_lung_edges_2_75_HumanNet_-1_1scaled_top_25000_edges_symbol.csv', offset=FALSE)
SWEET_SOW_lung <- calculate_SOW('./SWEET/SWEET_lung_no_doubles_no_loops_HN_-1_1scaled_top_25000_edges.csv', offset=FALSE)

sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,24)]
sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
# Remove sample 'SOW$`0775`' from metadata
sample_metadata <- sample_metadata[-which(sample_metadata$DepMap_ID == '0775'),]

plot_boxplot <- function(network, sample_metadata, filename){
  
  #network <- SWEET_SOW_lung
  
  network <- data.frame(t(network)); network$DepMap_ID <- rownames(network)
  # Merge the 'network' matrix and 'sample_metadata' data frame
  merged_data <- merge(sample_metadata, network, by.x = "DepMap_ID", by.y = 'DepMap_ID')
  
  # Reshape the data from wide to long format
  long_data <- gather(merged_data, feature, value, -(DepMap_ID:lineage_subtype))
  cbbPalette <- c("#DDCC77","#661100")
  # Create the boxplot with color mapping
  pdf(filename)
  
  print(ggplot(long_data, aes(x = DepMap_ID, y = value, fill = lineage_subtype)) +
          geom_boxplot(show.legend = FALSE) +
          scale_fill_manual(values = cbbPalette ) +
          labs(title = "",
               x = "Sample",
               y = "Sum of absolute edge weights") + 
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)))
    
  dev.off()
}

plot_boxplot(SWEET_SOW_lung, sample_metadata, '../results/Rebuttal/SWEET/SOW_boxplots/SWEET_lung.pdf')
plot_boxplot(SSN_SOW_lung, sample_metadata, '../results/Rebuttal/SWEET/SOW_boxplots/SSN_lung.pdf')
plot_boxplot(CSN_SOW_lung, sample_metadata, '../results/Rebuttal/SWEET/SOW_boxplots/CSN_lung.pdf')
plot_boxplot(iENA_SOW_lung, sample_metadata, '../results/Rebuttal/SWEET/SOW_boxplots/iENA_lung.pdf')
plot_boxplot(SSPGI_SOW_lung, sample_metadata, '../results/Rebuttal/SWEET/SOW_boxplots/SSPGI_lung.pdf')
plot_boxplot(LIONESS_SOW_lung, sample_metadata, '../results/Rebuttal/SWEET/SOW_boxplots/LIONESS_lung.pdf')


############################################################################################################################
#### Brain
############################################################################################################################

plot_boxplot <- function(network, sample_metadata, filename){ # Same function but use other colors for plot
  
  #network <- SWEET_SOW_lung
  
  network <- data.frame(t(network)); network$DepMap_ID <- rownames(network)
  # Merge the 'network' matrix and 'sample_metadata' data frame
  merged_data <- merge(sample_metadata, network, by.x = "DepMap_ID", by.y = 'DepMap_ID')
  
  # Reshape the data from wide to long format
  long_data <- gather(merged_data, feature, value, -(DepMap_ID:lineage_subtype))
  cbbPalette <- c("#E69F00", "#009E73")
  # Create the boxplot with color mapping
  pdf(filename)
  
  print(ggplot(long_data, aes(x = DepMap_ID, y = value, fill = lineage_subtype)) +
          geom_boxplot(show.legend = FALSE) +
          scale_fill_manual(values = cbbPalette ) +
          labs(title = "",
               x = "Sample",
               y = "Sum of absolute edge weights") + 
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)))
  
  dev.off()
}

LIONESS_SOW_brain <- calculate_SOW('./LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset=FALSE)
SSN_SOW_brain  <- calculate_SOW('./SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset=FALSE)
CSN_SOW_brain  <-calculate_SOW('./CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv', offset=FALSE)
iENA_SOW_brain  <-calculate_SOW('./ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv',offset=FALSE)
SSPGI_SOW_brain  <-calculate_SOW('./SSPGI/SSPGI_brain_edges_2_75_HumanNet_-1_1scaled_top_25000_edges_symbol.csv', offset=FALSE)
SWEET_SOW_brain <- calculate_SOW('./SWEET/SWEET_brain_no_doubles_no_loops_HN_-1_1scaled_top_25000_edges.csv', offset=FALSE)

sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,19)]
sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))

# Select glioblastoma and medulloblastoma samples
sample_metadata <- sample_metadata[sample_metadata$Subtype %in% c('Glioblastoma', 'Medulloblastoma'),]
colnames(sample_metadata)[2] <- 'lineage_subtype'

plot_boxplot(SWEET_SOW_brain, sample_metadata, '../results/Rebuttal/SWEET/SOW_boxplots/SWEET_brain.pdf')
plot_boxplot(SSN_SOW_brain, sample_metadata, '../results/Rebuttal/SWEET/SOW_boxplots/SSN_brain.pdf')
plot_boxplot(CSN_SOW_brain, sample_metadata, '../results/Rebuttal/SWEET/SOW_boxplots/CSN_brain.pdf')
plot_boxplot(iENA_SOW_brain, sample_metadata, '../results/Rebuttal/SWEET/SOW_boxplots/iENA_brain.pdf')
plot_boxplot(SSPGI_SOW_brain, sample_metadata, '../results/Rebuttal/SWEET/SOW_boxplots/SSPGI_brain.pdf')
plot_boxplot(LIONESS_SOW_brain, sample_metadata, '../results/Rebuttal/SWEET/SOW_boxplots/LIONESS_brain.pdf')

