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

setwd('/home/boris/Documents/PhD/single_sample/')
#############################################################################################################################
#### Proteomics data analysis
#############################################################################################################################

proteomics <- fread('./protein_quant_current_normalized.csv')
# sample_info <- fread('./Table_S1_Sample_Information.xlsx', fill=TRUE)
sample_metadata <- fread('./networks/20Q4_v2_sample_info.csv')
sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))

samples <- names(fread('./networks/LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', header=TRUE, data.table=FALSE, fill=TRUE))[-c(1,2)]
sample_metadata <- sample_metadata[sample_metadata$DepMap_ID %in% samples, ]
proteomics <- data.frame(proteomics[, c(2,49:426)])
colnames(proteomics) <- str_replace(colnames(proteomics), ('_Ten.*$'),'')
proteomics <- dplyr::select(proteomics,-c('CAL120_BREAST', 'SW948_LARGE_INTESTINE', 'HCT15_LARGE_INTESTINE')) # Remove non-unique cell lines
# 12 197 out of 12 755 gene symbols are unique, remove the remaining ones
proteomics <- data.frame(proteomics %>% group_by(Gene_Symbol) %>% filter(row_number()==1))
rownames(proteomics) <- proteomics$Gene_Symbol; proteomics$Gene_Symbol <- NULL
proteomics <- proteomics[, colnames(proteomics) %in% sample_metadata$CCLE_Name] #26 matching samples
colnames(proteomics) <- sample_metadata$DepMap_ID[match(colnames(proteomics), sample_metadata$CCLE_Name)]
sample_metadata <- data.frame(sample_metadata[sample_metadata$DepMap_ID %in% colnames(proteomics), c(1,24)])
# Check normalization
proteomics %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))

# What to do with missing values? First try just setting them to zero
#proteomics[is.na(proteomics)] <- 0 
proteomics <- proteomics[, order(names(proteomics))]

#############################################################################################################################
#### Transcriptomics data analysis
#############################################################################################################################

network <- fread('./networks/LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet.csv', header = TRUE, data.table = FALSE, fill = TRUE)
genes <- unique(as.vector(as.matrix(network[,c(1,2)])))

expr <- fread('selection_lung_expression_data.csv', header=TRUE, fill = TRUE, data.table = FALSE)
expr$V1 <- substr(expr$V1, 7, nchar(expr$V1))
rownames(expr) <- expr$V1; expr$V1 <- NULL
expr <- as.data.frame(t(expr))
# expr <- expr[rownames(expr) %in% genes, ] # This selects for genes present in the aggregate network (which are highly variable)... Maybe no filtering is better?
# rownames(expr) <- str_replace(row.names(expr), " \\(.*", "")
expr$gene <- str_replace(row.names(expr), " \\(.*", "")

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host='https://www.ensembl.org')
id_ensembl <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'gene_biotype'),
                    filters='hgnc_symbol', values=expr$gene, mart=ensembl)

RNAseq <- merge(expr, id_ensembl, by.x = 'gene', by.y = 'hgnc_symbol')
RNAseq_coding <- as.data.frame(RNAseq[!duplicated(RNAseq$ensembl_gene_id), ]) # 5 non-unique ensembl ID's
rownames(RNAseq_coding) <- RNAseq_coding$ensembl_gene_id
RNAseq_coding <- RNAseq_coding[RNAseq_coding$gene_biotype == 'protein_coding', ]
RNAseq_coding <- as.data.frame(RNAseq_coding %>% group_by(gene) %>% filter(row_number()==1)) # Just keep the first entry, they have the same counts anyway
rownames(RNAseq_coding) <- RNAseq_coding$gene
RNAseq_coding$gene_biotype <- NULL; RNAseq_coding$gene <- NULL; RNAseq_coding$ensembl_gene_id <- NULL


#############################################################################################################################
#### TRy gene-wise correlations across samples 
#############################################################################################################################

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

# LIONESS_SOW <- calculate_SOW('./LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset=TRUE)
# SSN_SOW <- calculate_SOW('./SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv', offset=TRUE)
# CSN_SOW <-calculate_SOW('./CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset=TRUE)
# iENA_SOW <-calculate_SOW('./ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv',offset=TRUE)
# SSPGI_SOW <-calculate_SOW('./SSPGI/SSPGI_lung_edges_2_75_HumanNet_top_25000_edges_symbol.csv', offset=TRUE)

LIONESS_SOW <- calculate_SOW('./LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset=FALSE)
SSN_SOW <- calculate_SOW('./SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', offset=FALSE)
CSN_SOW <-calculate_SOW('./CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv', offset=FALSE)
iENA_SOW <-calculate_SOW('./ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv',offset=FALSE)
SSPGI_SOW <-calculate_SOW('./SSPGI/SSPGI_lung_edges_2_75_HumanNet_-1_1scaled_top_25000_edges_symbol.csv', offset=FALSE)
SWEET_SOW <- calculate_SOW('./SWEET/SWEET_lung_no_doubles_no_loops_HN_-1_1scaled_top_25000_edges.csv', offset=FALSE)

proteomics <- na.omit(proteomics) # 5620 proteins without NA values in the 28 coupled samples
RNAseq_coding <- RNAseq_coding[rownames(RNAseq_coding) %in% rownames(proteomics), colnames(RNAseq_coding) %in% colnames(proteomics)]

gene_wise_corr <- function(proteomics_dat, expr_SOW, aggregate=FALSE){
  # expr_SOW <- LIONESS_SOW
  # proteomics_dat <- proteomics
  dat <- expr_SOW[rownames(expr_SOW) %in% rownames(proteomics_dat), colnames(expr_SOW) %in% colnames(proteomics_dat)]; proteomics_dat <- proteomics_dat[rownames(proteomics_dat) %in% rownames(expr_SOW), ]
  
  proteomics_dat <- proteomics_dat[order(rownames(proteomics_dat)), ]
  ord <- match(rownames(proteomics_dat), rownames(dat))
  dat <- dat[ord, ]
  
  dat <- data.frame(t(dat))
  proteomics_dat <- data.frame(t(proteomics_dat))
  
  # Empty dataframe to store correlation score per gene
  corr_scores <- data.frame(matrix(nrow=nrow(dat), ncol=1))
  rownames(corr_scores) <- rownames(proteomics_dat); colnames(corr_scores) <- 'corr_coeff'
  
  for (i in 1:nrow(corr_scores)){
    # i <- 3
    prot_dat <- as.vector(proteomics_dat[i,])
    expr_dat <- as.vector(dat[i,])
    corr_res <- cor.test(as.numeric(prot_dat), as.numeric(expr_dat), method='spearman')
    corr_scores[i,] <- corr_res$estimate
  }
  return(corr_scores)
}

expression <- gene_wise_corr(proteomics_dat = proteomics, expr_SOW = RNAseq_coding); expression$gene <- rownames(expression)
LIONESS <- gene_wise_corr(proteomics, LIONESS_SOW); LIONESS$gene <- rownames(LIONESS)
SSN <- gene_wise_corr(proteomics, SSN_SOW); SSN$gene <- rownames(SSN)
CSN <- gene_wise_corr(proteomics, CSN_SOW); CSN$gene <- rownames(CSN)
iENA <- gene_wise_corr(proteomics, iENA_SOW); iENA$gene <- rownames(iENA)
SSPGI <- gene_wise_corr(proteomics, SSPGI_SOW); SSPGI$gene <- rownames(SSPGI)
SWEET <- gene_wise_corr(proteomics, SWEET_SOW); SWEET$gene <- rownames(SWEET)

aggregate <-fread('./aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', data.table = FALSE, header=TRUE, fill=TRUE)
aggregate$reg <- str_replace(aggregate$reg,  "\\(.*", ""); aggregate$tar <- str_replace(aggregate$tar,  "\\(.*", "")
genes <- unique(as.vector(as.matrix(aggregate[,c(1,2)])))
aggregate_SOW <- data.frame(matrix(nrow = length(genes), ncol = 1))
rownames(aggregate_SOW) <- genes; colnames(aggregate_SOW) <- 'aggregate'

for (i in 1:length(genes)){
  #i <- 1
  gene <- rownames(aggregate_SOW)[i]
  ind <- which(aggregate[,1] %in% gene)
  ind <- c(ind, which(aggregate[,2] %in% gene))
  # network[,-c(1,2)] <- lapply(network[,-c(1,2)], )
  
  # SOW[i,] <- colSums(abs(network[ind, -c(1,2)])) # This is the approach used in the basic figure now included in the paper, comment lines 36:42, 55 out to change
  aggregate_SOW[i,] <- sum(network[ind, -c(1,2)])
}
# proteomics_mean <- rowMeans(proteomics); proteomics_mean <- proteomics_mean[names(proteomics_mean) %in% rownames(aggregate_SOW)]
proteomics <- proteomics[rownames(proteomics) %in% rownames(aggregate_SOW), ]; aggregate_SOW <- aggregate_SOW[rownames(aggregate_SOW) %in% rownames(proteomics), ]
aggregate_res <- vector('list', length=ncol(proteomics)); names(aggregate_res) <- colnames(proteomics)
for (i in 1:ncol(proteomics)){
  aggregate_res[i] <- cor.test(aggregate_SOW, proteomics[,i], method = 'spearman')$estimate
}
aggregate_res <- data.frame(t(data.frame(aggregate_res))); rownames(aggregate_res) <- gsub('X','', rownames(aggregate_res))
aggregate_res$gene <- rownames(aggregate_res); colnames(aggregate_res) <- c('corr_coeff', 'gene')
df_list <- list(expression, SSN, LIONESS, SWEET, iENA, CSN, SSPGI, aggregate_res)

merged <- data.frame(df_list %>% purrr::reduce(left_join, by='gene'))

rownames(merged) <- merged$gene; merged$gene <- NULL
colnames(merged) <- c('Expression', 'SSN', 'LIONESS', 'SWEET', 'iENA', 'CSN', 'SSPGI','Aggregate'); merged <- abs(merged)
merged_long <- reshape2::melt(merged)

kruskal.test(value~variable ,data = merged_long) 
res <- dunnTest(value~variable, data=merged_long, method = 'holm')
#write.table(res, file = '../results/Rebuttal/SWEET/proteomics_sample_wise_correlation_lung_stats.txt', sep = '\t', quote = FALSE, row.names = FALSE)

# Should I use ANOVA (parametric) or Kruskal-Wallis test (non-parametric) here?
cbbPalette <- c("#FFFFFF", "#a9a9a9", "#a9a9a9", "#a9a9a9", "#a9a9a9", "#a9a9a9", "#a9a9a9", "#000000")
comparisons <- list(#c('Expression', 'SSN'), c('Expression', 'LIONESS'), c('Expression', 'iENA'), c('Expression', 'CSN'), c('Expression', 'SSPGI'), c('Expression', 'Aggregate'),
  c('Aggregate', 'SSPGI'), c('Aggregate', 'CSN'), c('Aggregate', 'iENA'), c('Aggregate', 'SWEET'), c('Aggregate', 'LIONESS'), c('Aggregate', 'SSN'),
  c('iENA', 'LIONESS'), c('LIONESS', 'SSPGI'), c('SSN', 'SSPGI'))
# c('Aggregate', 'SSPGI'), c('Aggregate', 'CSN'),  c('Aggregate', 'iENA'), c('Aggregate', 'LIONESS'), c('Aggregate', 'SSN')
# c('Aggregate', 'SSN'), c('Aggregate', 'LIONESS'), c('Aggregate', 'iENA'), c('Aggregate', 'CSN'), c('Aggregate', 'SSPGI')
pdf('../results/Rebuttal/SWEET/proteomics_sample_wise_correlation_lung_scaled_withSWEET.pdf')
ggplot(merged_long, aes(x = variable, y=value, fill = variable)) + geom_boxplot() + 
  stat_compare_means(comparisons = comparisons, label = 'p.signif', vjust=0.75, label.x = 5, label.y = c(0.3, 0.32,0.34,0.36,0.38, 0.4, 0.42, 0.44, 0.46), tip.length = 0.02) +
  stat_compare_means(label.x = 5.75, label.y = 0.260) +
  ggtitle('') + scale_fill_manual(values=cbbPalette) +
  xlab("") + ylab("Correlation coefficient") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=0, hjust=0.5, size = 12),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        # panel.grid.minor = element_line(size=0.2, linetype='solid', colour='black'), 
        # panel.grid.major = element_line(size=0.1, linetype='solid', colour='grey'),
        axis.line = element_line(size = 0.5, linetype = 'solid', colour = 'black'),
        legend.position = 'none') 
  
dev.off()



# Figure for OncoPoint
df_list <- list(SSN, LIONESS, iENA, CSN, SSPGI, aggregate_res)
merged <- data.frame(df_list %>% reduce(left_join, by='gene'))
rownames(merged) <- merged$gene; merged$gene <- NULL
colnames(merged) <- c('SSN', 'LIONESS', 'iENA', 'CSN', 'SSPGI','Aggregate'); merged <- abs(merged)
merged_long <- reshape2::melt(merged)
comparisons <- list(c('Aggregate', 'SSPGI'), c('Aggregate', 'CSN'),  c('Aggregate', 'iENA'), c('Aggregate', 'LIONESS'), c('Aggregate', 'SSN'))
cbbPalette <- c("#a9a9a9", "#a9a9a9", "#a9a9a9", "#a9a9a9", "#a9a9a9", "#FFFFFF")

pdf('/home/boris/shared/proteomics_lung_noExpr_poster.pdf')
ggplot(merged_long, aes(x = variable, y=value, fill = variable)) + geom_boxplot() + 
  stat_compare_means(comparisons = comparisons, label = 'p.signif', vjust=0.75, label.x = 5, label.y = c(0.3, 0.32,0.34,0.36,0.38), tip.length = 0.02) +
  stat_compare_means(label.x = 4.75, label.y = 0.260) +
  ggtitle('Correlation between protein abundance and node importance \n in sample-specific networks') + scale_fill_manual(values=cbbPalette) +
  xlab("") + ylab("Mean correlation coefficient") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle=0, hjust=0.5, size = 12),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        # panel.grid.minor = element_line(size=0.2, linetype='solid', colour='black'), 
        # panel.grid.major = element_line(size=0.1, linetype='solid', colour='grey'),
        axis.line = element_line(size = 0.5, linetype = 'solid', colour = 'black'),
        legend.position = 'none') 
dev.off()


