#!/usr/bin/Rscript

##################################################################################################################
##### Limma analysis on expression data
##################################################################################################################

library(data.table)
library(stringr)
library(limma)
library(EnhancedVolcano)
library(plyr)
library(edgeR)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpubr)
library(jamba)
library(data.table)
library(ggplot2)
library(ggfortify)
library(ggpubr)
library(biomaRt)
setwd('/home/boris/Documents/PhD/single_sample/')
# I will select for genes that are also present in the networks -> making things more comparable
network <- fread('./networks/LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet.csv', header = TRUE, data.table = FALSE, fill = TRUE)
genes <- unique(as.vector(as.matrix(network[,c(1,2)])))

expr <- fread('selection_lung_expression_data.csv', header=TRUE, fill = TRUE, data.table = FALSE)
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


sample_metadata <- fread('./networks/20Q4_v2_sample_info.csv')[,c(1,24)]
sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
# Remove the carcinoid sample
RNAseq_coding$`0775` <- NULL
sample_metadata <- data.frame(sample_metadata[sample_metadata$DepMap_ID %in% colnames(RNAseq_coding)])
rownames(sample_metadata) <- sample_metadata$DepMap_ID; sample_metadata[,1] <- NULL

# Perform differential analysis on the entire dataset first, use limma with emperical Bayes (also used in Lopes-Ramos 2021) 
design <- as.formula(~ 0 + lineage_subtype)
design <- model.matrix(design, data = sample_metadata); colnames(design) <- c('NSCLC', 'SCLC')

# Filter out genes with low counts
expr <- DGEList(RNAseq_coding)
keep <- filterByExpr(expr, design)
expr <- expr[keep, , keep.lib.sizes=FALSE]
expr <- calcNormFactors(expr)
logcpm <- cpm(expr, log=TRUE, prior.count =1)

fit <- lmFit(logcpm, design)
cont.matrix <- makeContrasts(NSCLC-SCLC, levels = design) #This does NSCLC vs SCLC, so positive LFCs are higher in NSCLC compared to SCLC
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)
res <- topTable(fit, adjust = "BH", number = nrow(expr))
write.table(res, 'Lung_DGE_analysis.txt', quote = FALSE)

# Check one gene
gene <- 'SYT4'
scores <- data.frame(logcpm[gene, ])
scores <- merge(scores, sample_metadata, by=0)
rownames(scores) <- scores$Row.names; scores$Row.names <- NULL

# Create violing plot
sample_size = scores %>% group_by(lineage_subtype) %>% summarise(num=n()); colnames(sample_size) <- c('Lineage subtype', 'num')
colnames(scores) <- c('Count', 'Lineage subtype')

print(scores %>%
        left_join(sample_size) %>%
        mutate(myaxis = paste0(`Lineage subtype`, '\n', 'n=', num)) %>%
        ggplot(aes(x = myaxis, y = `Count`, fill=`Lineage subtype`)) +
        geom_violin(width=1.4) +
        geom_boxplot(width=0.1, color="grey", alpha=0.2) +
        scale_fill_viridis(discrete = TRUE) +
        theme_light() +
        theme(legend.position="none", plot.title = element_text(size=11, hjust =0.5)) +
        # ggtitle(paste0(out, '_', gene)) +
        xlab(""))

pdf(paste0('./results/DGE/lung/expression_lung_Volcano.pdf'))
print(EnhancedVolcano(res,
                      lab = rownames(res),
                      x = 'logFC',
                      y = 'adj.P.Val',
                      labSize = 3,
                      pCutoff = 0.05,
                      title = 'Lung expression data'))
dev.off()


############################################################################################################
#### GO analysis
############################################################################################################
setwd('/home/boris/Documents/PhD/single_sample/results/DGE/lung')
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(org.Hs.eg.db)

gene_list_all <- abs(res$logFC)
names(gene_list_all) <- row.names(res)
gene_list_all = sort(gene_list_all, decreasing = TRUE) # Sorting in decreasing order is required for clusterProfiler

# keyType This is the source of the annotation (gene ids). The options vary for each annotation. 
organism <- org.Hs.eg.db
keytypes(organism)
head(keys(organism, keytype="SYMBOL"))

gse_all <- gseGO(geneList=gene_list_all, 
                 ont ="ALL", # Which ontology: cellular component, molecular function, biological process
                 keyType = "SYMBOL", #genesymbol
                 # nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time) | The function itself suggests to not use the nPerm argument
                 minGSSize = 3, # Minimum number of genes for enrichment
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "none")
table <- gse_all@result[,c(1:10)]
fwrite(table, file='./expression_lung_gsea_all.txt')

# Create dotplot visualization of gsea analysis
pdf('genes_all_gsea_dotplot2.pdf')
dotplot(gse_all, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

gse_bio <- gseGO(geneList=gene_list_all, 
                 ont ="BP", # Which ontology: cellular component, molecular function, biological process
                 keyType = "SYMBOL", #genesymbol
                 nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                 minGSSize = 3, # Minimum number of genes for enrichment
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "none")

# Create dotplot visualization of gsea analysis
pdf('genes_BioProcess_gsea_dotplot.pdf')
dotplot(gse_bio, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

gse_mf <- gseGO(geneList=gene_list_all, 
                ont ="MF", # Which ontology: cellular component, molecular function, biological process
                keyType = "SYMBOL", #genesymbol
                nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                minGSSize = 3, # Minimum number of genes for enrichment
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = organism, 
                pAdjustMethod = "none")

# Create dotplot visualization of gsea analysis
pdf('genes_MolFunction_gsea_dotplot.pdf')
dotplot(gse_mf, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

# pdf('all_genes_enrichmentPerGOcategory.pdf')
# dotplot(gse_all, showCategory=10, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("Main") + theme(axis.text.y = element_text(size=6))
# dev.off()

gse_cc <- gseGO(geneList=gene_list_all, 
                ont ="CC", # Which ontology: cellular component, molecular function, biological process
                keyType = "SYMBOL", #genesymbol
                nPerm = 100, # Number of permutations, the more the more accurate the results (but takes time)
                minGSSize = 3, # Minimum number of genes for enrichment
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = organism, 
                pAdjustMethod = "none")

# Create dotplot visualization of gsea analysis
pdf('genes_CelComponent_gsea_dotplot.pdf')
dotplot(gse_cc, showCategory=10, split=".sign") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size=6))
dev.off()

# pdf('all_genes_enrichmentPerGOcategory.pdf')
# dotplot(gse_all, showCategory=10, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free") + ggtitle("Main") + theme(axis.text.y = element_text(size=6))
# dev.off()

termsim_all <- pairwise_termsim(gse_all)
pdf('./genes_all_enrichmentMap.pdf')
emapplot(termsim_all, showCategory = 10) #pie_scale=1.5,layout="kk"
dev.off()

edox <- setReadable(gse_all, 'org.Hs.eg.db', 'ENSEMBL')
pdf('./genes_all_cnetplot.pdf')
cnetplot(edox, categorySize="pvalue", foldChange=gene_list_all, showCategory = 10) #3
dev.off()

pdf('genes_all_enrichmentDistribution.pdf')
ridgeplot(gse_all) + labs(x = "enrichment distribution") + theme(axis.text.y = element_text(size=4))
dev.off()

pdf('./genes_all_Heatplot_down.pdf')
heatplot(gse_all, foldChange=gene_list_all, showCategory = 5) + theme(axis.text.x = element_text(size=2))
dev.off()

sign_genes <- res[(abs(res$logFC) >= 1 & res$adj.P.Val <= 0.05), ]
sign_up <- sign_genes[sign_genes$logFC >= 1, ]
sign_down <- sign_genes[sign_genes$logFC <= -1, ]
write.table(rownames(sign_up), './genes_up.txt', quote=FALSE, row.names = FALSE, sep = '\t')
write.table(rownames(sign_down), './genes_down.txt', quote=FALSE, row.names = FALSE, sep = '\t')


########################################################################################################################################################
#### PCA analysis
########################################################################################################################################################

plotPCA <- function(data, output, ranked = FALSE){
  # network <- './LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet.csv'
  # network <- fread(network, data.table = FALSE, header = TRUE, fill = TRUE)
  data <- data.frame(logcpm)
  colnames(data) <- gsub('X', '', colnames(data))
  output <- 'Expression data'
  # offset <- FALSE
  # ranked <- TRUE
  
  # Data needs to be transposed
  data <- as.data.frame(t(data))
  
  # Add sample annotations to this dataframe
  sample_metadata <- fread('./networks/20Q4_v2_sample_info.csv')[,c(1,24)]
  sample_metadata$DepMap_ID <- substr(sample_metadata$DepMap_ID, 7, nchar(sample_metadata$DepMap_ID))
  sample_metadata <- data.frame(sample_metadata[sample_metadata$DepMap_ID %in% rownames(data)])
  rownames(sample_metadata) <- sample_metadata$DepMap_ID; sample_metadata[,1] <- NULL
  
  data <- merge(data, sample_metadata, by=0)
  rownames(data) <- data$Row.names; data$Row.names <- NULL; data <- data[, c(ncol(data),2:ncol(data)-1)]
  
  # Set colorblind friendly colors
  # The palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  # To use for line and point colors, add
  scale_colour_manual(values=cbbPalette)
  
  
  PCA <- prcomp(data[, -c(1)], scale=FALSE)
  pca_plot <- autoplot(PCA, data = data, colour = 'lineage_subtype')
  pca_plot <- pca_plot + theme_light() + ggtitle(output) + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=cbbPalette)
  pdf('/home/boris/Documents/PhD/single_sample/results/DGE/lung/PCA.pdf')
  pca_plot
  dev.off()

  return(pca_plot)
}
PCA_plot <- plotPCA(data.frame(logcpm, 'Expression data'))

################################################################################
#### Try to create some kind of network visualization of the aggregate network
################################################################################
library(igraph)
aggregate_lung <- fread('/home/boris/Documents/PhD/single_sample/networks/aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv', header=TRUE, data.table=FALSE, fill=TRUE)
aggregate_lung$reg <- str_replace(aggregate_lung$reg, "\\(.*","")
aggregate_lung$tar <- str_replace(aggregate_lung$tar, "\\(.*","")


NSCLC_drivers <- fread("/home/boris/Documents/PhD/single_sample/networks/IntOGen-DriverGenes_NSCLC.tsv")
SCLC_drivers <- fread("/home/boris/Documents/PhD/single_sample/networks/IntOGen-DriverGenes_SCLC.tsv")

genes <- unique(union(NSCLC_drivers$Symbol, SCLC_drivers$Symbol))
aggregate_lung_filtered <- aggregate_lung[(aggregate_lung$reg %in% genes | aggregate_lung$tar %in% genes), ]

aggregate_lung_filtered$type <- 'PP'
aggregate_lung_filtered <- aggregate_lung_filtered[, c(1,4,2,3)]
write.table(aggregate_lung_filtered, '/home/boris/Documents/PhD/single_sample/networks/aggregate/lung_cytoscape.txt', sep = '\t', quote=FALSE, row.names = FALSE)

# colnames(aggregate_lung) <- c("From", "To", "Weight")
# network_graph <- graph.data.frame(aggregate_lung, directed = FALSE)
# degrees <- degree(network_graph)
# to_keep <- degrees[degrees >= 100]
# aggregate_lung_filtered <- aggregate_lung[(aggregate_lung$From %in% names(to_keep) | aggregate_lung$To %in% names(to_keep)), ]

# Create an attribute table that will allow to selectively color driver genes in Cytoscape networks

network_genes <- data.frame(unique(as.vector(as.matrix(aggregate_lung_filtered[,c(1,3)]))))
network_genes$node_type <- 'Gene'; colnames(network_genes) <- c('node_id', 'node_type')
network_genes$node_type[network_genes$node_id %in% genes] <-'Driver'
write.table(network_genes, '/home/boris/Documents/PhD/single_sample/networks/aggregate/lung_cytoscape_attributes.txt', sep = '\t', quote=FALSE, row.names = FALSE)



