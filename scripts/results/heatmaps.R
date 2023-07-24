#!/usr/bin/Rscript

library(data.table)
library(gtools)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(stringr)
library(DESeq2)
library(edgeR)
library(stringr)

setwd('/home/boris/Documents/PhD/single_sample/networks')


create_heatmap_brain <- function(input_df, outfile, network = FALSE, keep_index =FALSE){
  
  # input_df <- './SSPGI/SSPGI_brain_edges_2_75_HumanNet_top_25000_edges.csv'
  # network <- TRUE
  # outfile <- '/home/boris/Documents/PhD/single_sample/LIONESS_Heatmap_Z_score_Spearman_brain.jpeg'
  # keep_index <- FALSE
  
  expr_df <- fread(input_df, data.table=F, header=T)
  
  colnames(expr_df) <- str_replace(colnames(expr_df), "X", '')
  if (keep_index == FALSE){
    row.names(expr_df) <- expr_df[,1]
    expr_df[,1] <- NULL 
  }
  
  sample_names <- colnames(expr_df)
  expr_df
  
  if (network == TRUE){
    expr_df <- data.frame(expr_df[, -c(1,2)])
  }
  
  colnames(expr_df) <- str_replace(colnames(expr_df), "X", '')
  
  
  sample_info <- read.csv('/home/boris/Documents/PhD/single_sample/20Q4_v2_sample_info.csv', header=TRUE, sep=",", fill=TRUE)
  sample_info
  
  cell_lines <- paste("ACH-00", colnames(expr_df), sep="")
  cell_line_sample_info <- sample_info[match(cell_lines, sample_info$DepMap_ID) ,] 
  cell_line_sample_info
  
  spearman_matrix <- cor(expr_df, method="spearman")
  mu <- mean(spearman_matrix)
  sigma <- sd(spearman_matrix)
  zscore_matrix <- apply(spearman_matrix, c(1,2), function(x) (x-mu)/sigma)
  
  row.names(spearman_matrix) <-colnames(expr_df)
  colnames(spearman_matrix) <- colnames(expr_df)
  row.names(zscore_matrix) <-colnames(expr_df)
  colnames(zscore_matrix) <- colnames(expr_df)
  
  row_ha1 = rowAnnotation(disease_subtype = cell_line_sample_info$Subtype, show_annotation_name = FALSE, 
                          col = list(disease_subtype = c("Astrocytoma" = "#000000", "Glioblastoma" = "#E69F00", "Glioma" = "#56B4E9", 
                                                         "Medulloblastoma" = "#009E73", "Oligodendroglioma" = "#F0E442", "Meningioma"="#0072B2", 
                                                         "Primitive Neuroectodermal Tumor (PNET)"="#D55E00")), annotation_label="Disease subtype")
  
  
  col_fun = colorRamp2(c(-2, 0, 2), c("steelblue4", "white", "firebrick1")) # choose the colors for the heatmap 
  col_fun(seq(-16, 16)) #how much different colors should be visible 
  
  hm <- Heatmap(zscore_matrix,  name="Spearman correlation (Z-score)", 
                clustering_distance_rows = "spearman", clustering_distance_columns = "spearman", 
                clustering_method_rows = "ward.D",clustering_method_columns = "ward.D", width = unit(15, "cm"), 
                height = unit(15, "cm"), left_annotation = c(row_ha1), col=col_fun, row_dend_width = unit(2, "cm"), 
                row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7), show_column_dend = FALSE) # the plot 
  
  jpeg(paste0('../results/',outfile, '_Heatmap_Z_score_Spearman_brain.jpeg'), width=13, height=8, units="in", res=300)
  draw(hm, merge_legend = TRUE) # merge legend so all the legend parts are shown vertically  
  dev.off()
}

create_heatmap_brain('../JDS_brain_expression_2_75_var_only_names.csv', 'Expression')
create_heatmap_brain('./SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv', 'SSN', network = TRUE)
create_heatmap_brain('./LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv', 'LIONESS', network =TRUE)
create_heatmap_brain('./ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv', 'iENA', network =TRUE)
create_heatmap_brain('./CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv', 'CSN', network = TRUE, keep_index = TRUE)
create_heatmap_brain('./SSPGI/SSPGI_brain_edges_2_75_HumanNet_top_25000_edges_symbol.csv', 'SSPGI', network =TRUE)


create_heatmap_lung <- function(input_df, outfile, network = FALSE, keep_index =FALSE){
  
  # input_df <- './SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv'
  # network <- TRUE
  # outfile <- '/home/boris/Documents/PhD/single_sample/SSN_Heatmap_Z_score_Spearman_lung.jpeg'
  # keep_index <- FALSE
  
  expr_df <- fread(input_df, data.table=F, header=T)
  
  colnames(expr_df) <- str_replace(colnames(expr_df), "X", '')
  if (keep_index == FALSE){
    row.names(expr_df) <- expr_df[,1]
    expr_df[,1] <- NULL 
  }
  
  sample_names <- colnames(expr_df)
  expr_df
  
  if (network == TRUE){
    expr_df <- data.frame(expr_df[, -c(1,2)])
  }
  
  colnames(expr_df) <- str_replace(colnames(expr_df), "X", '')
  
  
  sample_info <- read.csv('/home/boris/Documents/PhD/single_sample/20Q4_v2_sample_info.csv', header=TRUE, sep=",", fill=TRUE)
  sample_info
  sample_info$lineage_subtype <- str_replace(sample_info$lineage_subtype, 'lung_carcinoid', 'Lung_carcinoid')
  
  cell_lines <- paste("ACH-00", colnames(expr_df), sep="")
  cell_line_sample_info <- sample_info[match(cell_lines, sample_info$DepMap_ID) ,] 
  cell_line_sample_info
  
  spearman_matrix <- cor(expr_df, method="spearman")
  mu <- mean(spearman_matrix)
  sigma <- sd(spearman_matrix)
  zscore_matrix <- apply(spearman_matrix, c(1,2), function(x) (x-mu)/sigma)
  
  row.names(spearman_matrix) <-colnames(expr_df)
  colnames(spearman_matrix) <- colnames(expr_df)
  row.names(zscore_matrix) <-colnames(expr_df)
  colnames(zscore_matrix) <- colnames(expr_df)
  
  row_ha1 = rowAnnotation(disease_subtype = cell_line_sample_info$lineage_subtype, show_annotation_name = FALSE, 
                          col = list(disease_subtype = c("Lung_carcinoid" = "#888888", "NSCLC" = "#DDCC77", "SCLC" = "#661100")),
                          annotation_label="Disease subtype")
  
  
  col_fun = colorRamp2(c(-2, 0, 2), c("steelblue4", "white", "firebrick1")) # choose the colors for the heatmap 
  col_fun(seq(-16, 16)) #how much different colors should be visible 
  
  hm <- Heatmap(zscore_matrix,  name="Spearman correlation (Z-score)", 
                clustering_distance_rows = "spearman", clustering_distance_columns = "spearman", 
                clustering_method_rows = "ward.D",clustering_method_columns = "ward.D", width = unit(15, "cm"), 
                height = unit(15, "cm"), left_annotation = c(row_ha1), col=col_fun, row_dend_width = unit(2, "cm"), 
                row_names_gp = gpar(fontsize = 7), column_names_gp = gpar(fontsize = 7), show_column_dend = FALSE) # the plot 
  
  jpeg(paste0('../results/',outfile, '_Heatmap_Z_score_Spearman_lung.jpeg'), width=13, height=8, units="in", res=300)
  draw(hm, merge_legend = TRUE) # merge legend so all the legend parts are shown vertically  
  dev.off()
}

create_heatmap_lung('../JDS_lung_expression_2_75_var_only_names.csv', 'Expression')
create_heatmap_lung('./SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv', 'SSN', network = TRUE)
create_heatmap_lung('./LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv', 'LIONESS', network =TRUE)
create_heatmap_lung('./ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv', 'iENA', network =TRUE)
create_heatmap_lung('./CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv', 'CSN', network = TRUE, keep_index = TRUE)
create_heatmap_lung('./SSPGI/SSPGI_lung_edges_2_75_HumanNet_top_25000_edges_symbol.csv', 'SSPGI', network =TRUE)

