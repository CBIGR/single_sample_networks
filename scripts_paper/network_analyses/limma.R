#!usr/bin/Rscript

library(data.table)
library(stringr)
library(limma)
library(EnhancedVolcano)
library(plyr)
library(edgeR)
library(ggplot2)
library(dplyr)
library(viridis)
library(jamba)
library(data.table)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(biomaRt)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
#library(hrbrthemes)
setwd('/home/boris/Documents/PhD/single_sample/networks')

Calculate_sum_of_edge_weights <- function(network, offset=FALSE, log=FALSE, abs=TRUE){
  
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

##############################################################################################################################
#### Perform limma analysis on sum of absolute edge weights (just like the Glass paper does on weights of individual edges)
##############################################################################################################################

Differential_nodes_lung <- function(networkfile, out, offset = FALSE, log=FALSE){
  
  # networkfile <- "./SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv"
  # offset = FALSE
  # out <- 'SSN'
  # log <- FALSE
  # abs <- FALSE
  
  network <- fread(networkfile, header=TRUE, data.table=FALSE, fill=TRUE)
  colnames(network) <- gsub('X', '', colnames(network))
  print(dim(network))
  if (offset == TRUE){ # Remove first column
    network[,1] <- NULL
  }
  
  # Extract table with sum of edge weights, rename them to gene symbols
  if (log == TRUE){ #SSPGI edge weights are way too large, results in a very high variance and no significant results
    SOW <- Calculate_sum_of_edge_weights(network, log = TRUE)
  } else {
    SOW <- Calculate_sum_of_edge_weights(network, abs=TRUE)
  }
  rownames(SOW) <- str_replace(row.names(SOW), "\\s.*", "")
  rownames(SOW) <- str_replace(row.names(SOW), "\\(.*", "")
  
  # Subset SOW dataframe to dataframe that only contains relevant edges, test for differences between 2 groups and create boxplots per driver gene
  sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,24)]
  sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  # Remove the carcinoid sample
  SOW$`0775` <- NULL
  sample_metadata <- data.frame(sample_metadata[sample_metadata$DepMap_ID %in% colnames(SOW)])
  rownames(sample_metadata) <- sample_metadata$DepMap_ID; sample_metadata[,1] <- NULL
  
  # Perform differential analysis on the entire dataset first, use limma with emperical Bayes (also used in Lopes-Ramos 2021) 
  design <- as.formula(~ 0 + lineage_subtype)
  design <- model.matrix(design, data = sample_metadata); colnames(design) <- c('NSCLC', 'SCLC')
  fit <- lmFit(SOW, design)
  cont.matrix <- makeContrasts(NSCLC-SCLC, levels = design)
  fit <- contrasts.fit(fit, cont.matrix)
  fit <- eBayes(fit)
  res <- topTable(fit, adjust = "BH", number = nrow(SOW))
  
  
  # # Plot results in a volcano plot
  pdf(paste0('../results/Rebuttal/SWEET/limma/', out, '_lung_Volcano_abs.pdf'))
  print(EnhancedVolcano(res,
                        lab = rownames(res),
                        x = 'logFC',
                        y = 'adj.P.Val',
                        labSize = 3,
                        title = out,
                        pCutoff = 0.05,
                        FCcutoff = 1))
  dev.off()

  return(res)
}

# networkfile <- './LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv'
# offset = FALSE
# log <- FALSE
# out <- 'LIONESS_scaled'

setwd('/home/boris/Documents/PhD/single_sample/networks')
Lioness_lung_diff_nodes <- Differential_nodes_lung('./LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', 'LIONESS', offset=FALSE)
SSN_lung_diff_nodes <- Differential_nodes_lung('./SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', 'SSN', offset=FALSE)
CSN_lung_diff_nodes <-Differential_nodes_lung('./CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv', 'CSN_scaled', offset=FALSE)
iENA_lung_diff_nodes <-Differential_nodes_lung('./ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', 'ssPCC', offset=FALSE)
SSPGI_lung_diff_nodes <-Differential_nodes_lung('./SSPGI/SSPGI_lung_edges_2_75_HumanNet_-1_1scaled_top_25000_edges_symbol.csv', 'SSPGI', offset=FALSE)
SWEET_lung_diff_nodes <-Differential_nodes_lung('./SWEET/SWEET_lung_no_doubles_no_loops_HN_-1_1scaled_top_25000_edges.csv', 'SWEET', offset=FALSE)


Differential_nodes_brain <- function(networkfile, out, offset = FALSE, log=FALSE){
  
  # networkfile <- './LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv'
  # log <- FALSE
  # offset <- FALSE
  
  network <- fread(networkfile, header=TRUE, fill=TRUE, data.table=FALSE)
  colnames(network) <- gsub('X', '', colnames(network))
  if (offset == TRUE){ # Remove first column
    network[,1] <- NULL
  }
  
  # Extract table with sum of edge weights, rename them to gene symbols
  if (log == TRUE){ #SSPGI edge weights are way too large, results in a very high variance and no significant results
    SOW <- Calculate_sum_of_edge_weights(network, log = TRUE)
  } else {
    SOW <- Calculate_sum_of_edge_weights(network, abs=TRUE)
  }
  rownames(SOW) <- str_replace(row.names(SOW), "\\s.*", "")
  rownames(SOW) <- str_replace(row.names(SOW), "\\(.*", "")
  
  # Subset SOW dataframe to dataframe that only contains relevant edges, test for differences between 2 groups and create boxplots per driver gene
  sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,19)]
  sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  sample_metadata <- data.frame(sample_metadata[sample_metadata$DepMap_ID %in% colnames(SOW)])
  # Select medulloblastoma and glioblastoma samples (the two largest groups)
  sample_metadata <- sample_metadata[sample_metadata$Subtype %in% c('Glioblastoma', 'Medulloblastoma'), ]
  rownames(sample_metadata) <- sample_metadata$DepMap_ID; sample_metadata$DepMap_ID <- NULL
  SOW <- SOW[, colnames(SOW) %in% row.names(sample_metadata)]
  
  # Perform differential analysis on the entire dataset first, use limma with emperical Bayes (also used in Lopes-Ramos 2021) 
  design <- as.formula(~ 0 + Subtype)
  design <- model.matrix(design, data = sample_metadata); colnames(design) <- c('Glioblastoma', 'Medulloblastoma')
  fit <- lmFit(SOW, design)
  cont.matrix <- makeContrasts(Glioblastoma-Medulloblastoma, levels = design)
  fit <- contrasts.fit(fit, cont.matrix)
  fit <- eBayes(fit)
  res <- topTable(fit, adjust = "BH", number = nrow(SOW))
  
  # Plot results in a volcano plot
  pdf(paste0('../results/Rebuttal/SWEET/limma/', out, '_brain_Volcano_abs.pdf'))
  print(EnhancedVolcano(res,
                        lab = rownames(res),
                        x = 'logFC',
                        y = 'adj.P.Val',
                        labSize = 3,
                        title = out,
                        pCutoff = 0.05,
                        FCcutoff = 1))
  dev.off()
  
 
  return(res)
}

setwd('/home/boris/Documents/PhD/single_sample/networks')
Lioness_brain_diff_nodes <- Differential_nodes_brain('./LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', 'LIONESS', offset=FALSE)
SSN_brain_diff_nodes <- Differential_nodes_brain('./SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', 'SSN', offset=FALSE)
CSN_brain_diff_nodes <-Differential_nodes_brain('./CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv', 'CSN_scaled', offset=FALSE)
iENA_brain_diff_nodes <-Differential_nodes_brain('./ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', 'ssPCC', offset=FALSE)
SSPGI_brain_diff_nodes <-Differential_nodes_brain('./SSPGI/SSPGI_brain_edges_2_75_HumanNet_-1_1scaled_top_25000_edges_symbol.csv', 'SSPGI', offset=FALSE)
SWEET_brain_diff_nodes <- Differential_nodes_brain('./SWEET/SWEET_brain_no_doubles_no_loops_HN_-1_1scaled_top_25000_edges.csv', 'SWEET', offset=FALSE)



####################################################################################################################################################################################
#### Test enrichment for IntOGen/COSMIC drivers
####################################################################################################################################################################################

SCLC_drivers <- fread("/home/boris/Documents/PhD/single_sample/networks/IntOGen-DriverGenes_SCLC.tsv")[,1]
NSCLC_drivers <- fread("/home/boris/Documents/PhD/single_sample/networks/IntOGen-DriverGenes_NSCLC.tsv")[,1]
lung_drivers <- rbind(SCLC_drivers, NSCLC_drivers); lung_drivers <- unique(as.vector(lung_drivers$Symbol))
GBM_drivers <- fread("/home/boris/Documents/PhD/single_sample/networks/IntOGen-DriverGenes_GBM.tsv")[,1]
MBL_drivers <- fread("/home/boris/Documents/PhD/single_sample/networks/IntOGen-DriverGenes_MBL.tsv")[,1]
brain_drivers <- rbind(MBL_drivers, MBL_drivers); brain_drivers <- unique(as.vector(brain_drivers$Symbol))
total_lung <- 5454
total_brain <- 4741

test_enrichment <- function(res, drivers, total, out){
  # res <- Lioness_lung_diff_nodes
  sign_res <- res[(abs(res$logFC) >= 1 & res$adj.P.Val <= 0.05), ]
  gene_set <- row.names(sign_res); print(length(gene_set))
  # drivers <- lung_drivers
  # total <- total_lung
  # out <- 'LIONESS'
  succes <- intersect(gene_set, drivers)
  possibilities <- length(drivers)
  total_genes <- total
  number_of_diff_nodes <- length(gene_set)
  
  test <- phyper(length(succes), possibilities, total_genes, number_of_diff_nodes)
  # pdf(paste0('../results/scaled/limma/', out, '_brain_Volcano_Abs.pdf'))
  # print(EnhancedVolcano(res,
  #                       lab = rownames(res),
  #                       x = 'logFC',
  #                       y = 'adj.P.Val',
  #                       labSize = 3,
  #                       title = out,
  #                       subtitle = test,
  #                       pCutoff = 0.05,
  #                       FCcutoff = 1))
  # dev.off()
  return(test)
}

test_enrichment(Lioness_lung_diff_nodes, lung_drivers, total_lung, out='LIONESS')
test_enrichment(SSN_lung_diff_nodes, lung_drivers, total_lung, out='SSN')
test_enrichment(CSN_lung_diff_nodes, lung_drivers, total_lung, out='CSN')
test_enrichment(iENA_lung_diff_nodes, lung_drivers, total_lung, out='iENA')
test_enrichment(SSPGI_lung_diff_nodes, lung_drivers, total_lung, out='SSPGI')
test_enrichment(SWEET_lung_diff_nodes, lung_drivers, total_lung, out='SWEET')

test_enrichment(Lioness_brain_diff_nodes, brain_drivers, total_brain, out='LIONESS')
test_enrichment(SSN_brain_diff_nodes, brain_drivers, total_brain, out='SSN')
test_enrichment(CSN_brain_diff_nodes, brain_drivers, total_brain, out='CSN')
test_enrichment(iENA_brain_diff_nodes, brain_drivers, total_brain, out='iENA')
test_enrichment(SSPGI_brain_diff_nodes, brain_drivers, total_brain, out='SSPGI')
test_enrichment(SWEET_brain_diff_nodes, brain_drivers, total_brain, out='SWEET')



####################################################################################################################################################################################
#### Test enrichment for Cell Model Passport drivers
####################################################################################################################################################################################

# Read in sample metadata and group samples based on subtype
sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,24,25)]
sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
# SCLC samples have no annotation in sub_subtype, just put the subtype there
sample_metadata$lineage_sub_subtype[sample_metadata$lineage_sub_subtype == ""] <- sample_metadata$lineage_subtype[sample_metadata$lineage_sub_subtype == ""]
base <- '/home/boris/Documents/PhD/single_sample/networks/'
samples_in_networks <- colnames(fread(paste0(base, "LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), header = TRUE))[4:89]
sample_metadata <- sample_metadata[sample_metadata$DepMap_ID %in% samples_in_networks, ]

SCLC_samples <- sample_metadata[sample_metadata$lineage_sub_subtype == 'SCLC', ]
NSCLC_samples <- sample_metadata[sample_metadata$lineage_subtype == 'NSCLC', ]
MBL_samples <- sample_metadata[sample_metadata$Subtype == "Medulloblastoma"]
GBM_samples <- sample_metadata[sample_metadata$Subtype == "Glioblastoma"]

# Read in driver mutation data
driver_mutations <- fread('../drivers/Cellpassport_mutations.csv', header =TRUE, fill =TRUE, data.table = FALSE)
driver_mutations$DepMap_ID <- substr(driver_mutations$DepMap_ID, 7, nchar(driver_mutations$DepMap_ID))
colnames(driver_mutations) <- c('DepMap_ID', 'Subtype', 'lineage_subtype', 'lineage_sub_subtype', 'model_id', 'Symbol')


SCLC_drivers <- driver_mutations[driver_mutations$DepMap_ID %in% SCLC_samples$DepMap_ID, ]
NSCLC_drivers <- driver_mutations[driver_mutations$DepMap_ID %in% NSCLC_samples$DepMap_ID, ]
lung_drivers <- rbind(SCLC_drivers, NSCLC_drivers); lung_drivers <- unique(as.vector(lung_drivers$Symbol))
MBL_drivers <- driver_mutations[driver_mutations$DepMap_ID %in% MBL_samples$DepMap_ID, ]
GBM_drivers <- driver_mutations[driver_mutations$DepMap_ID %in% GBM_samples$DepMap_ID, ]
brain_drivers <- rbind(MBL_drivers, MBL_drivers); brain_drivers <- unique(as.vector(brain_drivers$Symbol))

total_lung <- 5454
total_brain <- 4741


test_enrichment(Lioness_lung_diff_nodes, lung_drivers, total_lung, out='LIONESS')
test_enrichment(SSN_lung_diff_nodes, lung_drivers, total_lung, out='SSN')
test_enrichment(CSN_lung_diff_nodes, lung_drivers, total_lung, out='CSN')
test_enrichment(iENA_lung_diff_nodes, lung_drivers, total_lung, out='iENA')
test_enrichment(SSPGI_lung_diff_nodes, lung_drivers, total_lung, out='SSPGI')
test_enrichment(SWEET_lung_diff_nodes, lung_drivers, total_lung, out='SWEET')

test_enrichment(Lioness_brain_diff_nodes, brain_drivers, total_brain, out='LIONESS')
test_enrichment(SSN_brain_diff_nodes, brain_drivers, total_brain, out='SSN')
test_enrichment(CSN_brain_diff_nodes, brain_drivers, total_brain, out='CSN')
test_enrichment(iENA_brain_diff_nodes, brain_drivers, total_brain, out='iENA')
test_enrichment(SSPGI_brain_diff_nodes, brain_drivers, total_brain, out='SSPGI')
test_enrichment(SWEET_brain_diff_nodes, brain_drivers, total_brain, out='SWEET')
