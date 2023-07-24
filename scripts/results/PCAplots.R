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
  network <- network[, -c(1,2)] 
  if (ranked == TRUE){
    network[] <- lapply(-network, rank, ties.method="min") 
  }
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
  PCA <- prcomp(network[, -c(1)], scale=FALSE)
  pca_plot <- autoplot(PCA, data = network, colour = 'Subtype')
  pca_plot <- pca_plot + theme_light() + ggtitle(output) + theme(plot.title = element_text(hjust = 0.5)) + 
    scale_colour_manual(values=cbbPalette) + xlim(c(-0.3,1)) + ylim(c(-0.6,0.5))
  pca_plot
  
  
  results <- PCA
  #calculate total variance explained by each principal component
  var_explained = results$sdev^2 / sum(results$sdev^2)
  
  #create scree plot
  library(ggplot2)
  
  pdf(paste0('/home/boris/Documents/PhD/single_sample/results/scaled/PCA_scree_lung', output, '.pdf'))
  print(qplot(c(1:86), var_explained) + 
    geom_line() + 
    xlab("Principal Component") + 
    ylab("Variance Explained") +
    ggtitle(output) +
    ylim(0, 1))
  dev.off()
  
  results <- PCA$x
  # fviz_nbclust(results, FUNcluster=kmeans, k.max = 8)
  # fviz_nbclust(results, FUNcluster=kmeans, method="gap_stat", k.max = 8)+ theme_classic()
  
  # k-means clustering on PCA results
  km1<-eclust(results, "kmeans", hc_metric="euclidean",k=3)
  
  sample_metadata[sample_metadata$Subtype == 'NSCLC', ] <- 1
  sample_metadata[sample_metadata$Subtype == 'SCLC', ] <- 2
  sample_metadata[sample_metadata$Subtype == 'Lung_carcinoid', ] <- 3
  
  kclusters <- km1$cluster
  correct_clusters <- sample_metadata$Subtype
  
  library(fossil)
  adj.rand <- adj.rand.index(unname(kclusters), correct_clusters); print(adj.rand)
  
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

plotPCA_lung_expr <- function(expr, offset = FALSE, output, ranked = FALSE){
  # network <- './LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv'
  # output <- 'LIONESS'
  # offset <- FALSE
  # ranked <- FALSE
  
  #expr <- '/home/boris/Documents/PhD/single_sample/selection_lung_expression_data.csv'
  
  expr <- fread(expr, header=TRUE, fill = TRUE, data.table = FALSE)
  expr$V1 <- substr(expr$V1, 7, nchar(expr$V1))
  rownames(expr) <- expr$V1; expr$V1 <- NULL
  expr <- as.data.frame(t(expr))
  # expr <- expr[rownames(expr) %in% genes, ] # This selects for genes present in the aggregate network (which are highly variable)... Maybe no filtering is better?
  # rownames(expr) <- str_replace(row.names(expr), " \\(.*", "")
  expr$gene <- str_replace(row.names(expr), " \\(.*", "")
  
  
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host='www.ensembl.org')
  
  id_ensembl <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'gene_biotype'),
                      filters='hgnc_symbol', values=expr$gene, mart=ensembl)
  
  RNAseq <- merge(expr, id_ensembl, by.x = 'gene', by.y = 'hgnc_symbol')
  RNAseq_coding <- as.data.frame(RNAseq[!duplicated(RNAseq$ensembl_gene_id), ]) # 5 non-unique ensembl ID's
  rownames(RNAseq_coding) <- RNAseq_coding$ensembl_gene_id
  RNAseq_coding <- RNAseq_coding[RNAseq_coding$gene_biotype == 'protein_coding', ]
  RNAseq_coding <- as.data.frame(RNAseq_coding %>% group_by(gene) %>% filter(row_number()==1)) # Just keep the first entry, they have the same counts anyway
  rownames(RNAseq_coding) <- RNAseq_coding$gene
  RNAseq_coding$gene_biotype <- NULL; RNAseq_coding$gene <- NULL; RNAseq_coding$ensembl_gene_id <- NULL
  
  # Add sample annotations to this dataframe
  sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,24)]
  sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  sample_metadata <- data.frame(sample_metadata[sample_metadata$DepMap_ID %in% colnames(expr)])
  rownames(sample_metadata) <- sample_metadata$DepMap_ID; sample_metadata[,1] <- NULL
  colnames(sample_metadata) <- c('Subtype'); sample_metadata$Subtype <- str_replace(sample_metadata$Subtype, 'lung_carcinoid', 'Lung_carcinoid')
  expr <- as.data.frame(t(RNAseq_coding))
  
  expr <- merge(expr, sample_metadata, by=0)
  rownames(expr) <- expr$Row.names; expr$Row.names <- NULL; expr <- expr[, c(ncol(expr),2:ncol(expr)-1)]
  #cbbPalette <- c("#009E73", "#D55E00", "#000000") # Colorblind palette
  cbbPalette <- c("#888888", "#DDCC77","#661100")
  # Set columns to numeric
  expr[-c(1)] <- sapply(expr[-c(1)],as.numeric)
  
  PCA <- prcomp(expr[, -c(1)], scale=FALSE)
  pca_plot <- autoplot(PCA, data = expr, colour = 'Subtype')
  pca_plot <- pca_plot + theme_light() + ggtitle('Expression') + theme(plot.title = element_text(hjust = 0.5)) + 
    scale_colour_manual(values=cbbPalette) + xlim(c(-0.3,1)) + ylim(c(-0.6,0.5))
  pca_plot
  
  results <- PCA
  #calculate total variance explained by each principal component
  var_explained = results$sdev^2 / sum(results$sdev^2)
  
  #create scree plot
  library(ggplot2)
  
  pdf(paste0('/home/boris/Documents/PhD/single_sample/results/scaled/PCA_scree_lung', output, '.pdf'))
  print(qplot(c(1:86), var_explained) + 
          geom_line() + 
          xlab("Principal Component") + 
          ylab("Variance Explained") +
          ggtitle(output) +
          ylim(0, 1))
  dev.off()
  
  results <- PCA$x
  # fviz_nbclust(results, FUNcluster=kmeans, k.max = 8)
  # fviz_nbclust(results, FUNcluster=kmeans, method="gap_stat", k.max = 8)+ theme_classic()
  
  # k-means clustering on PCA results
  km1<-eclust(results, "kmeans", hc_metric="euclidean",k=3)
  
  sample_metadata[sample_metadata$Subtype == 'NSCLC', ] <- 3
  sample_metadata[sample_metadata$Subtype == 'SCLC', ] <- 2
  sample_metadata[sample_metadata$Subtype == 'Lung_carcinoid', ] <- 1
  
  kclusters <- km1$cluster
  correct_clusters <- sample_metadata$Subtype
  
  library(fossil)
  adj.rand <- adj.rand.index(unname(kclusters), correct_clusters); print(adj.rand)
  
  return(pca_plot)
}

expression_lung <- plotPCA_lung_expr('/home/boris/Documents/PhD/single_sample/selection_lung_expression_data.csv', offset = FALSE, 'Expression')

leg <- get_legend(lioness_lung_scaled_top25k)
lioness_lung_scaled_top25k <- lioness_lung_scaled_top25k + theme(legend.position = "none")
SSN_lung_scaled_top25k <- SSN_lung_scaled_top25k + theme(legend.position = "none")
ssPCC_lung_scaled_top25k <- ssPCC_lung_scaled_top25k + theme(legend.position = "none")
SSPGI_lung_scaled_top25k <- SSPGI_lung_scaled_top25k + theme(legend.position = "none")
CSN_lung_scaled_top25k <- CSN_lung_scaled_top25k + theme(legend.position = "none")
plots <- list(SSN_lung_scaled_top25k, lioness_lung_scaled_top25k, ssPCC_lung_scaled_top25k, CSN_lung_scaled_top25k, SSPGI_lung_scaled_top25k, leg)
pdf('/home/boris/Documents/PhD/single_sample/results/scaled/PCA/PCA_top25k_lung_sameAxes.pdf')
print(ggarrange(plotlist = plots,
                labels = c("a", "b", "c", "d", "e"),
                ncol = 2, nrow = 3
))
dev.off()

expression_lung <- expression_lung + theme(legend.position = "none")
plots <- list(expression_lung, SSN_lung_scaled_top25k, lioness_lung_scaled_top25k, ssPCC_lung_scaled_top25k, CSN_lung_scaled_top25k, SSPGI_lung_scaled_top25k)
pdf('/home/boris/Documents/PhD/single_sample/results/scaled/PCA/PCA_top25k_lung_sameAxes_withExpr.pdf')
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
  network <- network[, -c(1,2)] 
  if (ranked == TRUE){
    network[] <- lapply(-network, rank, ties.method="min") 
  }
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
  PCA <- prcomp(network[, -c(1)], scale=FALSE)
  pca_plot <- autoplot(PCA, data = network, colour = 'Subtype')
  pca_plot <- pca_plot + theme_light() + ggtitle(output) + theme(plot.title = element_text(hjust = 0.5)) + 
    scale_colour_manual(values=cbbPalette) + xlim(c(-0.75, 1)) + ylim(c(-0.75,1))
  pca_plot
  
  results <- PCA
  #calculate total variance explained by each principal component
  var_explained = results$sdev^2 / sum(results$sdev^2)
  
  #create scree plot
  library(ggplot2)
  
  pdf(paste0('/home/boris/Documents/PhD/single_sample/results/scaled/PCA_scree_brain', output, '.pdf'))
  print(qplot(c(1:67), var_explained) + 
          geom_line() + 
          xlab("Principal Component") + 
          ylab("Variance Explained") +
          ggtitle(output) +
          ylim(0, 1))
  dev.off()
  
  sample_metadata$name <- rownames(sample_metadata)
  # Focus on glioblastoma and medulloblastoma samples
  to_keep <- rownames(sample_metadata[(sample_metadata$Subtype == 'Glioblastoma' | sample_metadata$Subtype == 'Medulloblastoma'), ])
  sample_metadata <- sample_metadata[to_keep, ]; sample_metadata$name <- NULL
  results <- data.frame(PCA$x[to_keep, ])
  # fviz_nbclust(results, FUNcluster=kmeans, k.max = 8)
  # fviz_nbclust(results, FUNcluster=kmeans, method="gap_stat", k.max = 8)+ theme_classic()
  # 
  # k-means clustering on PCA results
  km1<-eclust(results[,c(1,2)], "kmeans", hc_metric="euclidean",k=2)
  
  sample_metadata[sample_metadata$Subtype == 'Glioblastoma', ] <- 2
  sample_metadata[sample_metadata$Subtype == 'Medulloblastoma', ] <- 1
  
  kclusters <- km1$cluster
  correct_clusters <- sample_metadata$Subtype
  
  library(fossil)
  adj.rand <- adj.rand.index(unname(kclusters), correct_clusters); print(adj.rand)
  return(pca_plot)
}


# lioness_brain_scaled_top25k <- plotPCA_brain('./LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset = FALSE, 'LIONESS')
# CSN_brain_scaled_top25k <- plotPCA_brain('./CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset = FALSE, 'CSN')
# SSN_brain_scaled_top25k <- plotPCA_brain('./SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset = FALSE, 'SSN')
# ssPCC_brain_scaled_top25k <- plotPCA_brain('./ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset = FALSE, 'iENA')
# SSPGI_brain_scaled_top25k <- plotPCA_brain('./SSPGI/SSPGI_brain_edges_2_75_HumanNet_-1_1scaled_top_25000_edges_symbol.csv', offset = FALSE, 'SSPGI')

lioness_brain_scaled_top25k <- plotPCA_brain('./LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset = TRUE, 'LIONESS')
CSN_brain_scaled_top25k <- plotPCA_brain('./CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv', offset = FALSE, 'CSN')
SSN_brain_scaled_top25k <- plotPCA_brain('./SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset = TRUE, 'SSN')
ssPCC_brain_scaled_top25k <- plotPCA_brain('./ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset = TRUE, 'iENA')
SSPGI_brain_scaled_top25k <- plotPCA_brain('./SSPGI/SSPGI_brain_edges_2_75_HumanNet_top_25000_edges_symbol.csv', offset = TRUE, 'SSPGI')


plotPCA_brain_expr <- function(expr, offset = FALSE, output, ranked = FALSE){
  # network <- './LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv'
  # output <- 'LIONESS'
  # offset <- FALSE
  # ranked <- TRUE
  
  #expr <- '/home/boris/Documents/PhD/single_sample/selection_brain_expression_data.csv'
  
  expr <- fread(expr, header=TRUE, fill = TRUE, data.table = FALSE)
  expr$V1 <- substr(expr$V1, 7, nchar(expr$V1))
  rownames(expr) <- expr$V1; expr$V1 <- NULL
  expr <- as.data.frame(t(expr))
  # expr <- expr[rownames(expr) %in% genes, ] # This selects for genes present in the aggregate network (which are highly variable)... Maybe no filtering is better?
  # rownames(expr) <- str_replace(row.names(expr), " \\(.*", "")
  expr$gene <- str_replace(row.names(expr), " \\(.*", "")
  
  
  ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host='www.ensembl.org')
  
  id_ensembl <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'gene_biotype'),
                     filters='hgnc_symbol', values=expr$gene, mart=ensembl)
  
  RNAseq <- merge(expr, id_ensembl, by.x = 'gene', by.y = 'hgnc_symbol')
  RNAseq_coding <- as.data.frame(RNAseq[!duplicated(RNAseq$ensembl_gene_id), ]) # 5 non-unique ensembl ID's
  rownames(RNAseq_coding) <- RNAseq_coding$ensembl_gene_id
  RNAseq_coding <- RNAseq_coding[RNAseq_coding$gene_biotype == 'protein_coding', ]
  RNAseq_coding <- as.data.frame(RNAseq_coding %>% group_by(gene) %>% filter(row_number()==1)) # Just keep the first entry, they have the same counts anyway
  rownames(RNAseq_coding) <- RNAseq_coding$gene
  RNAseq_coding$gene_biotype <- NULL; RNAseq_coding$gene <- NULL; RNAseq_coding$ensembl_gene_id <- NULL
  
  expr <- RNAseq_coding
  
  # Add sample annotations to this dataframe
  sample_metadata <- fread('20Q4_v2_sample_info.csv')[,c(1,19)]
  sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  sample_metadata <- data.frame(sample_metadata[sample_metadata$DepMap_ID %in% colnames(expr)])
  rownames(sample_metadata) <- sample_metadata$DepMap_ID; sample_metadata[,1] <- NULL
  colnames(sample_metadata) <- c('Subtype')
  
  expr <- as.data.frame(t(expr))
  
  network <- merge(expr, sample_metadata, by=0)
  rownames(network) <- network$Row.names; network$Row.names <- NULL; network <- network[, c(ncol(network),2:ncol(network)-1)]
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # Colorblind palette
  PCA <- prcomp(network[, -c(1)], scale=FALSE)
  pca_plot <- autoplot(PCA, data = network, colour = 'Subtype')
  pca_plot <- pca_plot + theme_light() + ggtitle(output) + theme(plot.title = element_text(hjust = 0.5)) + 
    scale_colour_manual(values=cbbPalette) + xlim(c(-0.75, 1)) + ylim(c(-0.75,1))
  pca_plot
  
  results <- PCA
  #calculate total variance explained by each principal component
  var_explained = results$sdev^2 / sum(results$sdev^2)
  
  #create scree plot
  library(ggplot2)
  
  pdf(paste0('/home/boris/Documents/PhD/single_sample/results/scaled/PCA_scree_brain', output, '.pdf'))
  print(qplot(c(1:67), var_explained) + 
          geom_line() + 
          xlab("Principal Component") + 
          ylab("Variance Explained") +
          ggtitle(output) +
          ylim(0, 1))
  dev.off()
  
  sample_metadata$name <- rownames(sample_metadata)
  # Focus on glioblastoma and medulloblastoma samples
  to_keep <- rownames(sample_metadata[(sample_metadata$Subtype == 'Glioblastoma' | sample_metadata$Subtype == 'Medulloblastoma'), ])
  sample_metadata <- sample_metadata[to_keep, ]; sample_metadata$name <- NULL
  results <- data.frame(PCA$x[to_keep, ])
  # fviz_nbclust(results, FUNcluster=kmeans, k.max = 8)
  # fviz_nbclust(results, FUNcluster=kmeans, method="gap_stat", k.max = 8)+ theme_classic()
  # 
  # k-means clustering on PCA results
  km1<-eclust(results[,c(1,2)], "kmeans", hc_metric="euclidean",k=2)
  
  sample_metadata[sample_metadata$Subtype == 'Glioblastoma', ] <- 2
  sample_metadata[sample_metadata$Subtype == 'Medulloblastoma', ] <- 1
  
  kclusters <- km1$cluster
  correct_clusters <- sample_metadata$Subtype
  
  library(fossil)
  adj.rand <- adj.rand.index(unname(kclusters), correct_clusters); print(adj.rand)
  return(pca_plot)
}

expression_brain <- plotPCA_brain_expr('/home/boris/Documents/PhD/single_sample/selection_brain_expression_data.csv', offset = FALSE, 'Expression')

leg <- get_legend(lioness_brain_scaled_top25k)
lioness_brain_scaled_top25k <- lioness_brain_scaled_top25k + theme(legend.position = "none")
SSN_brain_scaled_top25k <- SSN_brain_scaled_top25k + theme(legend.position = "none")
ssPCC_brain_scaled_top25k <- ssPCC_brain_scaled_top25k + theme(legend.position = "none")
SSPGI_brain_scaled_top25k <- SSPGI_brain_scaled_top25k + theme(legend.position = "none")
CSN_brain_scaled_top25k <- CSN_brain_scaled_top25k + theme(legend.position = "none")
plots <- list(SSN_brain_scaled_top25k, lioness_brain_scaled_top25k, ssPCC_brain_scaled_top25k, CSN_brain_scaled_top25k, SSPGI_brain_scaled_top25k, leg)
pdf('/home/boris/Documents/PhD/single_sample/results/scaled/PCA/PCA_top25k_brain_sameAxes.pdf')
print(ggarrange(plotlist = plots,
                labels = c("a", "b", "c", "d", "e"),
                ncol = 2, nrow = 3
))
dev.off()

expression_brain <- expression_brain + theme(legend.position='none')
plots <- list(expression_brain, SSN_brain_scaled_top25k, lioness_brain_scaled_top25k, ssPCC_brain_scaled_top25k, CSN_brain_scaled_top25k, SSPGI_brain_scaled_top25k)
pdf('/home/boris/Documents/PhD/single_sample/results/scaled/PCA/PCA_top25k_brain_sameAxes_withExpr.pdf')
print(ggarrange(plotlist = plots,
                labels = c("a", "b", "c", "d", "e","f"),
                ncol = 2, nrow = 3
))
dev.off()

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

