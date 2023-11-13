#!/usr/bin/Rscript

############################################################################################################################
#### Goal: PCA on (ranked) HumanNet networks -> can we see a subtype clustering in certain methods?
############################################################################################################################

library(data.table)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(stringr)
library(factoextra)
library(cluster)
library(biomaRt)
library(dplyr)


setwd('/home/boris/Documents/PhD/single_sample/networks/')

Calculate_sum_of_edge_weights <- function(network, offset=FALSE, log=FALSE, abs=TRUE){
  #abs <- TRUE
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
    # network[,-c(1,2)] <- lapply(network[,-c(1,2)], )
    if (abs == TRUE){
      SOW[i,] <- colSums(abs(network[ind, -c(1,2)])) # This is the approach used in the basic figure now included in the paper, comment lines 36:42, 55 out to change
    } else {
      SOW[i,] <- colSums(network[ind, -c(1,2)])
    }
  }
  if (log == TRUE){
    SOW <- log(abs(SOW+1))
  }
  return(SOW)
}

plotPCA_lung <- function(network, offset = FALSE, output, ranked = FALSE){
  # network <- './LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv'
  # output <- 'LIONESS'
  # offset <- FALSE
  # ranked <- FALSE
  
  network <- fread(network, data.table = FALSE, header = TRUE, fill = TRUE)
  colnames(network) <- gsub('X', '', colnames(network))
  
  if (offset == TRUE){
    network <- network[, -1]
  }
  
  # Network needs to be transposed, also remove regulator and target column
  #network <- network[, -c(1,2)] 
  if (ranked == TRUE){
    network[] <- lapply(-network, rank, ties.method="min") 
  }
  
  network <- Calculate_sum_of_edge_weights(network, offset=FALSE, log=FALSE, abs=TRUE)
  network <- as.data.frame(t(network))
  
  # Add sample annotations to this dataframe
  sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,24)]
  sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  sample_metadata <- data.frame(sample_metadata[sample_metadata$DepMap_ID %in% rownames(network)])
  rownames(sample_metadata) <- sample_metadata$DepMap_ID; sample_metadata[,1] <- NULL
  colnames(sample_metadata) <- c('Subtype'); sample_metadata$Subtype <- str_replace(sample_metadata$Subtype, 'lung_carcinoid', 'Lung_carcinoid')
  
  
  network <- merge(network, sample_metadata, by=0)
  rownames(network) <- network$Row.names; network$Row.names <- NULL; network <- network[, c(ncol(network),2:ncol(network)-1)]
  #cbbPalette <- c("#009E73", "#D55E00", "#000000") # Colorblind palette
  cbbPalette <- c("#888888", "#DDCC77","#661100")
  if (output == 'CSN'){
    PCA <- prcomp(network[, -c(1)], scale=FALSE) # Cannot scale binary network
  } else {
    PCA <- prcomp(network[, -c(1)], scale=TRUE)
  }
  
  pca_plot <- autoplot(PCA, data = network, colour = 'Subtype')
  pca_plot <- pca_plot + theme_light() + ggtitle(output) + theme(plot.title = element_text(hjust = 0.5)) + 
    scale_colour_manual(values=cbbPalette) + xlim(c(-0.3,0.3)) + ylim(c(-0.4,0.5))
  pca_plot
  
  return(pca_plot)
}

#### Analysis on scaled top25k networks

# lioness_lung_scaled_top25k <- plotPCA_lung('./LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset = FALSE, 'LIONESS')
# CSN_lung_scaled_top25k <- plotPCA_lung('./CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset = FALSE, 'CSN')
# SSN_lung_scaled_top25k <- plotPCA_lung('./SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset = FALSE, 'SSN')
# ssPCC_lung_scaled_top25k <- plotPCA_lung('./ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset = FALSE, 'iENA')
# SSPGI_lung_scaled_top25k <- plotPCA_lung('./SSPGI/SSPGI_lung_edges_2_75_HumanNet_-1_1scaled_top_25000_edges_symbol.csv', offset = FALSE, 'SSPGI')

lioness_lung_scaled_top25k <- plotPCA_lung('./LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset = TRUE, 'LIONESS')
CSN_lung_scaled_top25k <- plotPCA_lung('./CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'CSN')
SSN_lung_scaled_top25k <- plotPCA_lung('./SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset = TRUE, 'SSN')
ssPCC_lung_scaled_top25k <- plotPCA_lung('./ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset = TRUE, 'iENA')
SSPGI_lung_scaled_top25k <- plotPCA_lung('./SSPGI/SSPGI_lung_edges_2_75_HumanNet_top_25000_edges_symbol.csv', offset = TRUE, 'SSPGI')
SWEET_lung_scaled_top25k <- plotPCA_lung('./SWEET/SWEET_lung_no_doubles_no_loops_HN_-1_1scaled_top_25000_edges.csv', offset = FALSE, 'SWEET')

leg <- get_legend(lioness_lung_scaled_top25k)
lioness_lung_scaled_top25k <- lioness_lung_scaled_top25k + theme(legend.position = "none")
SSN_lung_scaled_top25k <- SSN_lung_scaled_top25k + theme(legend.position = "none")
ssPCC_lung_scaled_top25k <- ssPCC_lung_scaled_top25k + theme(legend.position = "none")
SSPGI_lung_scaled_top25k <- SSPGI_lung_scaled_top25k + theme(legend.position = "none")
CSN_lung_scaled_top25k <- CSN_lung_scaled_top25k + theme(legend.position = "none")
SWEET_lung_scaled_top25k <- SWEET_lung_scaled_top25k + theme(legend.position = "none")

plots <- list(SSN_lung_scaled_top25k, lioness_lung_scaled_top25k,SWEET_lung_scaled_top25k, ssPCC_lung_scaled_top25k, CSN_lung_scaled_top25k, SSPGI_lung_scaled_top25k)
pdf('/home/boris/Documents/PhD/single_sample/results/Rebuttal/PCA/PCA_top25k_lung_sameAxes_SOW_scaled_withSWEET.pdf')
print(ggarrange(plotlist = plots,
                labels = c("a", "b", "c", "d", "e", "f"),
                ncol = 2, nrow = 3
))
dev.off()

plotPCA_lung_subsubtypes <- function(network, offset = FALSE, output, ranked = FALSE){
  # network <- './LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv'
  # output <- 'LIONESS'
  # offset <- FALSE
  # ranked <- FALSE

  network <- fread(network, data.table = FALSE, header = TRUE, fill = TRUE)
  colnames(network) <- gsub('X', '', colnames(network))
  
  if (offset == TRUE){
    network <- network[, -1]
  }
  
  # Network needs to be transposed, also remove regulator and target column
  #network <- network[, -c(1,2)] 
  if (ranked == TRUE){
    network[] <- lapply(-network, rank, ties.method="min") 
  }
  
  network <- Calculate_sum_of_edge_weights(network, offset=FALSE, log=FALSE, abs=TRUE)
  network <- as.data.frame(t(network))
  
  # Add sample annotations to this dataframe
  sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,24,25)]
  #sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,19)]
  sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  sample_metadata <- data.frame(sample_metadata[sample_metadata$DepMap_ID %in% rownames(network)])
  sample_metadata$lineage_sub_subtype[sample_metadata$lineage_sub_subtype == ""] <- sample_metadata$lineage_subtype[sample_metadata$lineage_sub_subtype == ""]
  rownames(sample_metadata) <- sample_metadata$DepMap_ID; sample_metadata[,1] <- NULL; sample_metadata[,1] <- NULL
  colnames(sample_metadata) <- c('Subtype')
  #sample_metadata$Subtype <- str_replace(sample_metadata$Subtype, 'lung_carcinoid', 'Lung_carcinoid')
  
  sample_metadata$name <- rownames(sample_metadata)
  sample_metadata$name <- NULL
  
  
  network <- merge(network, sample_metadata, by=0)
  rownames(network) <- network$Row.names; network$Row.names <- NULL; network <- network[, c(ncol(network),2:ncol(network)-1)]
  #cbbPalette <- c("#009E73", "#D55E00", "#000000") # Colorblind palette
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  if (output == 'CSN'){
    PCA <- prcomp(network[, -c(1)], scale=FALSE) # Cannot scale binary network
  } else {
    PCA <- prcomp(network[, -c(1)], scale=TRUE)
  }
  
  pca_plot <- autoplot(PCA, data = network, colour = 'Subtype')
  pca_plot <- pca_plot + theme_light() + ggtitle(output) + theme(plot.title = element_text(hjust = 0.5)) + 
    scale_colour_manual(values=cbbPalette) + xlim(c(-0.3,0.3)) + ylim(c(-0.4,0.5))
  pca_plot
  
  return(pca_plot)
}

lioness_lung_scaled_top25k <- plotPCA_lung_subsubtypes('./LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset = TRUE, 'LIONESS')
CSN_lung_scaled_top25k <- plotPCA_lung_subsubtypes('./CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'CSN')
SSN_lung_scaled_top25k <- plotPCA_lung_subsubtypes('./SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset = TRUE, 'SSN')
ssPCC_lung_scaled_top25k <- plotPCA_lung_subsubtypes('./ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset = TRUE, 'iENA')
SSPGI_lung_scaled_top25k <- plotPCA_lung_subsubtypes('./SSPGI/SSPGI_lung_edges_2_75_HumanNet_top_25000_edges_symbol.csv', offset = TRUE, 'SSPGI')
SWEET_lung_scaled_top25k <- plotPCA_lung_subsubtypes('./SWEET/SWEET_lung_no_doubles_no_loops_HN_-1_1scaled_top_25000_edges.csv', offset = FALSE, 'SWEET')

leg <- get_legend(lioness_lung_scaled_top25k)
lioness_lung_scaled_top25k <- lioness_lung_scaled_top25k + theme(legend.position = "none")
SSN_lung_scaled_top25k <- SSN_lung_scaled_top25k + theme(legend.position = "none")
ssPCC_lung_scaled_top25k <- ssPCC_lung_scaled_top25k + theme(legend.position = "none")
SSPGI_lung_scaled_top25k <- SSPGI_lung_scaled_top25k + theme(legend.position = "none")
CSN_lung_scaled_top25k <- CSN_lung_scaled_top25k + theme(legend.position = "none")
SWEET_lung_scaled_top25k <- SWEET_lung_scaled_top25k + theme(legend.position = "none")

plots <- list(SSN_lung_scaled_top25k, lioness_lung_scaled_top25k,SWEET_lung_scaled_top25k, ssPCC_lung_scaled_top25k, CSN_lung_scaled_top25k, SSPGI_lung_scaled_top25k)
pdf('/home/boris/Documents/PhD/single_sample/results/Rebuttal/SWEET/PCA/PCA_top25k_lung_sameAxes_SOW_scaled_subsubtypes_withSWEET.pdf')
print(ggarrange(plotlist = plots,
                labels = c("a", "b", "c", "d", "e", "f"),
                ncol = 2, nrow = 3
))
dev.off()



plotPCA_brain <- function(network, offset = FALSE, output, ranked = FALSE){
  # network <- './LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv'
  # output <- 'LIONESS'
  # offset <- FALSE
  # ranked <- TRUE
  
  network <- fread(network, data.table = FALSE, header = TRUE, fill = TRUE)
  colnames(network) <- gsub('X', '', colnames(network))
  
  if (offset == TRUE){
    network <- network[, -1]
  }
  
  # Network needs to be transposed, also remove regulator and target column
  #network <- network[, -c(1,2)] 
  if (ranked == TRUE){
    network[] <- lapply(-network, rank, ties.method="min") 
  }
  
  network <- Calculate_sum_of_edge_weights(network, offset=FALSE, log=FALSE, abs=TRUE)
  network <- as.data.frame(t(network))
  
  # Add sample annotations to this dataframe
  sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,19)]
  sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  sample_metadata <- data.frame(sample_metadata[sample_metadata$DepMap_ID %in% rownames(network)])
  rownames(sample_metadata) <- sample_metadata$DepMap_ID; sample_metadata[,1] <- NULL
  colnames(sample_metadata) <- c('Subtype')
  
  network <- merge(network, sample_metadata, by=0)
  rownames(network) <- network$Row.names; network$Row.names <- NULL; network <- network[, c(ncol(network),2:ncol(network)-1)]
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # Colorblind palette
  if (output == 'CSN'){
    PCA <- prcomp(network[, -c(1)], scale=FALSE) # Cannot scale binary network
  } else {
    PCA <- prcomp(network[, -c(1)], scale=TRUE)
  }
  pca_plot <- autoplot(PCA, data = network, colour = 'Subtype')
  pca_plot <- pca_plot + theme_light() + ggtitle(output) + theme(plot.title = element_text(hjust = 0.5)) + 
    scale_colour_manual(values=cbbPalette) + xlim(c(-0.75, 1)) + ylim(c(-0.75,1))
  pca_plot
  
  return(pca_plot)
}



lioness_brain_scaled_top25k <- plotPCA_brain('./LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset = TRUE, 'LIONESS')
CSN_brain_scaled_top25k <- plotPCA_brain('./CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'CSN')
SSN_brain_scaled_top25k <- plotPCA_brain('./SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset = TRUE, 'SSN')
ssPCC_brain_scaled_top25k <- plotPCA_brain('./ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset = TRUE, 'iENA')
SSPGI_brain_scaled_top25k <- plotPCA_brain('./SSPGI/SSPGI_brain_edges_2_75_HumanNet_top_25000_edges_symbol.csv', offset = TRUE, 'SSPGI')
SWEET_brain_scaled_top25k <- plotPCA_brain('./SWEET/SWEET_brain_no_doubles_no_loops_HN_-1_1scaled_top_25000_edges.csv', offset = FALSE, 'SWEET')



leg <- get_legend(lioness_brain_scaled_top25k)
lioness_brain_scaled_top25k <- lioness_brain_scaled_top25k + theme(legend.position = "none")
SSN_brain_scaled_top25k <- SSN_brain_scaled_top25k + theme(legend.position = "none")
ssPCC_brain_scaled_top25k <- ssPCC_brain_scaled_top25k + theme(legend.position = "none")
SSPGI_brain_scaled_top25k <- SSPGI_brain_scaled_top25k + theme(legend.position = "none")
CSN_brain_scaled_top25k <- CSN_brain_scaled_top25k + theme(legend.position = "none")
SWEET_brain_scaled_top25k <- SWEET_brain_scaled_top25k + theme(legend.position = "none")
plots <- list(SSN_brain_scaled_top25k, lioness_brain_scaled_top25k, ssPCC_brain_scaled_top25k, CSN_brain_scaled_top25k, SSPGI_brain_scaled_top25k, leg)
pdf('/home/boris/Documents/PhD/single_sample/results/Rebuttal/PCA/PCA_top25k_brain_sameAxes_SOW_scaled.pdf')
print(ggarrange(plotlist = plots,
                labels = c("a", "b", "c", "d", "e"),
                ncol = 2, nrow = 3
))
dev.off()

plots <- list(SSN_brain_scaled_top25k, lioness_brain_scaled_top25k, SWEET_brain_scaled_top25k, ssPCC_brain_scaled_top25k, CSN_brain_scaled_top25k, SSPGI_brain_scaled_top25k)
pdf('/home/boris/Documents/PhD/single_sample/results/Rebuttal/PCA/PCA_top25k_brain_sameAxes_SOW_scaled_withSWEET.pdf')
print(ggarrange(plotlist = plots,
                labels = c("a", "b", "c", "d", "e"),
                ncol = 2, nrow = 3
))
dev.off()

# expression_brain <- expression_brain + theme(legend.position='none')
# plots <- list(expression_brain, SSN_brain_scaled_top25k, lioness_brain_scaled_top25k, ssPCC_brain_scaled_top25k, CSN_brain_scaled_top25k, SSPGI_brain_scaled_top25k)
# pdf('/home/boris/Documents/PhD/single_sample/results/scaled/PCA/PCA_top25k_brain_sameAxes_withExpr.pdf')
# print(ggarrange(plotlist = plots,
#                 labels = c("a", "b", "c", "d", "e","f"),
#                 ncol = 2, nrow = 3
# ))
# dev.off()

# #### Analysis on normal HumanNet networks (no scaling, no top nodes)
# 
# lioness_lung_scaled_top25k <- plotPCA('./LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'LIONESS')
# CSN_lung_scaled_top25k <- plotPCA('./CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'CSN')
# SSN_lung_scaled_top25k <- plotPCA('./SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'SSN')
# ssPCC_lung_scaled_top25k <- plotPCA('./ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'iENA')
# SSPGI_lung_scaled_top25k <- plotPCA('./SSPGI/SSPGI_lung_edges_2_75_HumanNet.csv', offset = TRUE, 'SSPGI')
# 
# leg <- get_legend(lioness_lung_scaled_top25k)
# lioness_lung_scaled_top25k <- lioness_lung_scaled_top25k + theme(legend.position = "none")
# SSN_lung_scaled_top25k <- SSN_lung_scaled_top25k + theme(legend.position = "none")
# ssPCC_lung_scaled_top25k <- ssPCC_lung_scaled_top25k + theme(legend.position = "none")
# SSPGI_lung_scaled_top25k <- SSPGI_lung_scaled_top25k + theme(legend.position = "none")
# CSN_lung_scaled_top25k <- CSN_lung_scaled_top25k + theme(legend.position = 'none')
# plots <- list(lioness_lung_scaled_top25k, SSN_lung_scaled_top25k, ssPCC_lung_scaled_top25k, SSPGI_lung_scaled_top25k, CSN_lung_scaled_top25k, leg)
# pdf('/home/boris/Documents/PhD/single_sample/results/PCA/PCA_HumanNet_lung.pdf')
# print(ggarrange(plotlist = plots,
#                 labels = c("A", "B", "C", "D", "E"),
#                 ncol = 2, nrow = 3
# ))
# dev.off()
# 
# 
# lioness_brain_scaled_top25k <- plotPCA('./LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'LIONESS')
# CSN_brain_scaled_top25k <- plotPCA('./CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'CSN')
# SSN_brain_scaled_top25k <- plotPCA('./SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'SSN')
# ssPCC_brain_scaled_top25k <- plotPCA('./ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'iENA')
# SSPGI_brain_scaled_top25k <- plotPCA('./SSPGI/SSPGI_brain_edges_2_75_HumanNet.csv', offset = TRUE, 'SSPGI')
# 
# leg <- get_legend(lioness_brain_scaled_top25k)
# lioness_brain_scaled_top25k <- lioness_brain_scaled_top25k + theme(legend.position = "none")
# SSN_brain_scaled_top25k <- SSN_brain_scaled_top25k + theme(legend.position = "none")
# ssPCC_brain_scaled_top25k <- ssPCC_brain_scaled_top25k + theme(legend.position = "none")
# SSPGI_brain_scaled_top25k <- SSPGI_brain_scaled_top25k + theme(legend.position = "none")
# CSN_brain_scaled_top25k <- CSN_brain_scaled_top25k + theme(legend.position = 'none')
# plots <- list(lioness_brain_scaled_top25k, SSN_brain_scaled_top25k, ssPCC_brain_scaled_top25k, SSPGI_brain_scaled_top25k, CSN_brain_scaled_top25k,  leg)
# pdf('/home/boris/Documents/PhD/single_sample/results/PCA/PCA_HumanNet_brain.pdf')
# print(ggarrange(plotlist = plots,
#                 labels = c("A", "B", "C", "D", "E"),
#                 ncol = 2, nrow = 3
# ))
# dev.off()
# 
# 
# #### Analysis on ranked HUmanNet networks
# lioness_lung_scaled_top25k <- plotPCA('./LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'LIONESS', ranked =TRUE)
# CSN_lung_scaled_top25k <- plotPCA('./CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'CSN')
# SSN_lung_scaled_top25k <- plotPCA('./SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'SSN', ranked =TRUE)
# ssPCC_lung_scaled_top25k <- plotPCA('./ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'iENA', ranked =TRUE)
# SSPGI_lung_scaled_top25k <- plotPCA('./SSPGI/SSPGI_lung_edges_2_75_HumanNet.csv', offset = TRUE, 'SSPGI', ranked =TRUE)
# 
# leg <- get_legend(lioness_lung_scaled_top25k)
# lioness_lung_scaled_top25k <- lioness_lung_scaled_top25k + theme(legend.position = "none")
# SSN_lung_scaled_top25k <- SSN_lung_scaled_top25k + theme(legend.position = "none")
# ssPCC_lung_scaled_top25k <- ssPCC_lung_scaled_top25k + theme(legend.position = "none")
# SSPGI_lung_scaled_top25k <- SSPGI_lung_scaled_top25k + theme(legend.position = "none")
# CSN_lung_scaled_top25k <- CSN_lung_scaled_top25k + theme(legend.position = 'none')
# plots <- list(lioness_lung_scaled_top25k, SSN_lung_scaled_top25k, ssPCC_lung_scaled_top25k, SSPGI_lung_scaled_top25k, CSN_lung_scaled_top25k, leg)
# pdf('/home/boris/Documents/PhD/single_sample/results/PCA/PCA_HumanNet_ranked_lung.pdf')
# print(ggarrange(plotlist = plots,
#                 labels = c("A", "B", "C", "D", "E"),
#                 ncol = 2, nrow = 3
# ))
# dev.off()
# 
# 
# lioness_brain_scaled_top25k <- plotPCA('./LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'LIONESS', ranked =TRUE)
# CSN_brain_scaled_top25k <- plotPCA('./CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'CSN')
# SSN_brain_scaled_top25k <- plotPCA('./SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'SSN', ranked =TRUE)
# ssPCC_brain_scaled_top25k <- plotPCA('./ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'iENA', ranked =TRUE)
# SSPGI_brain_scaled_top25k <- plotPCA('./SSPGI/SSPGI_brain_edges_2_75_HumanNet.csv', offset = TRUE, 'SSPGI', ranked =TRUE)
# 
# leg <- get_legend(lioness_brain_scaled_top25k)
# lioness_brain_scaled_top25k <- lioness_brain_scaled_top25k + theme(legend.position = "none")
# SSN_brain_scaled_top25k <- SSN_brain_scaled_top25k + theme(legend.position = "none")
# ssPCC_brain_scaled_top25k <- ssPCC_brain_scaled_top25k + theme(legend.position = "none")
# SSPGI_brain_scaled_top25k <- SSPGI_brain_scaled_top25k + theme(legend.position = "none")
# CSN_brain_scaled_top25k <- CSN_brain_scaled_top25k + theme(legend.position = 'none')
# plots <- list(lioness_brain_scaled_top25k, SSN_brain_scaled_top25k, ssPCC_brain_scaled_top25k, SSPGI_brain_scaled_top25k, CSN_brain_scaled_top25k,  leg)
# pdf('/home/boris/Documents/PhD/single_sample/results/PCA/PCA_HumanNet_ranked_brain.pdf')
# print(ggarrange(plotlist = plots,
#                 labels = c("A", "B", "C", "D", "E"),
#                 ncol = 2, nrow = 3
# ))
# dev.off()

