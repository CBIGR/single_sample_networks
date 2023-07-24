
# save.image('../results/scaled/limma/differential_nodes_brain_all.RData')
# 
# # Get nodes with differential SOW for Enrichr analysis
# Lioness_brain_diff_nodes <- Lioness_brain_diff_nodes[(Lioness_brain_diff_nodes$adj.P.Val <= 0.05 & abs(Lioness_brain_diff_nodes$logFC) >= 1), ]
# write.table(row.names(Lioness_brain_diff_nodes), file = '../results/scaled/limma/Lioness_brain_scaled_diff_SOWnodes.txt', row.names = FALSE, quote = FALSE)
# 
# SSN_brain_diff_nodes <- SSN_brain_diff_nodes[(SSN_brain_diff_nodes$adj.P.Val <= 0.05 & abs(SSN_brain_diff_nodes$logFC) >= 1), ]
# write.table(row.names(SSN_brain_diff_nodes), file = '../results/scaled/limma/SSN_brain_scaled_diff_SOWnodes.txt', row.names = FALSE, quote = FALSE)
# 
# CSN_brain_diff_nodes <- CSN_brain_diff_nodes[(CSN_brain_diff_nodes$adj.P.Val <= 0.05 & abs(CSN_brain_diff_nodes$logFC) >= 1), ]
# write.table(row.names(CSN_brain_diff_nodes), file = '../results/scaled/limma/CSN_brain_scaled_diff_SOWnodes.txt', row.names = FALSE, quote = FALSE)
# 
# iENA_brain_diff_nodes <- iENA_brain_diff_nodes[(iENA_brain_diff_nodes$adj.P.Val <= 0.05 & abs(iENA_brain_diff_nodes$logFC) >= 1), ]
# write.table(row.names(iENA_brain_diff_nodes), file = '../results/scaled/limma/iENA_brain_scaled_diff_SOWnodes.txt', row.names = FALSE, quote = FALSE)
# 
# SSPGI_brain_diff_nodes <- SSPGI_brain_diff_nodes[(SSPGI_brain_diff_nodes$adj.P.Val <= 0.05 & abs(SSPGI_brain_diff_nodes$logFC) >= 1), ]
# write.table(row.names(SSPGI_brain_diff_nodes), file = '../results/scaled/limma/SSPGI_brain_scaled_diff_SOWnodes.txt', row.names = FALSE, quote = FALSE)
# 


####################################################################################################################################################################################
#### Is there a significant enrichment of subtype specific driver genes vs all known genes (or vs general cancer driver genes or vs the total number of genes present in a network)
####################################################################################################################################################################################

setwd('/home/boris/Documents/PhD/single_sample/networks')

# Now check enrichment for known drivers
SCLC_drivers_intogen <- fread('./IntOGen-DriverGenes_SCLC.tsv')[,1]
NSCLC_drivers_intogen <- fread('./IntOGen-DriverGenes_NSCLC.tsv')[,1]
SCLC_drivers_census <- fread('./Census_SCLC_19_06_2023.tsv')[,1]
NSCLC_drivers_census <- fread('./Census_NSCLC_19_06_2023.tsv')[,1]

SCLC_drivers <- rbind(SCLC_drivers_intogen, SCLC_drivers_census, use.names=FALSE)
NSCLC_drivers <- rbind(NSCLC_drivers_intogen, NSCLC_drivers_census, use.names=FALSE)
lung_drivers <- rbind(SCLC_drivers, NSCLC_drivers); lung_drivers <- unique(as.vector(lung_drivers$Symbol))


# Now check enrichment for known drivers
MBL_drivers_intogen <- fread('./IntOGen-DriverGenes_MBL.tsv')[,1]
HGG_drivers_intogen <- fread('./IntOGen-DriverGenes_GBM.tsv')[,1]
MBL_drivers_census <- fread('./Census_MBL_19_06_2023.tsv')[,1]
HGG_drivers_census <- fread('./Census_GBM_19_06_2023.tsv')[,1]

MBL_drivers <- rbind(MBL_drivers_intogen, MBL_drivers_census, use.names=FALSE)
HGG_drivers <- rbind(HGG_drivers_intogen, HGG_drivers_census, use.names=FALSE)
brain_drivers <- rbind(HGG_drivers, MBL_drivers); brain_drivers <- unique(as.vector(brain_drivers$Symbol))

convertToSymbol <- function(genelist){
  genelist <- str_replace(genelist, "\\s.*","")
  genelist <- str_replace(genelist, "\\(.*", "")
  return(genelist)
}

test_enrichment <- function(res, drivers, total, out, lung=TRUE){
  
  # res <- Lioness_lung_diff_nodes
  # drivers <- lung_drivers
  # # total <- total_lung
  # out <- 'LIONESS'
  # lung <- TRUE
  
  sign_res <- res[(abs(res$logFC) >= 1 & res$adj.P.Val <= 0.05), ]
  gene_set <- row.names(sign_res)
  
  if (lung == TRUE){
    aggregate <- fread('./aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv', fill = TRUE, header = TRUE, data.table = FALSE) #Should we use the complete aggregate network or the top 25k edges?
  } else {
    aggregate <- fread('./aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet.csv', fill = TRUE, header = TRUE, data.table = FALSE) #Should we use the complete aggregate network or the top 25k edges?
  }
  
  aggregate$reg <- convertToSymbol(aggregate$reg)
  aggregate$tar <- convertToSymbol(aggregate$tar)
  HumanNet_genes <- unique(as.vector(as.matrix(aggregate[,c(1,2)])))
  
  #succes <- intersect(gene_set, drivers)
  drivers_in_network <- intersect(drivers, HumanNet_genes)
  number_of_diff_nodes <- length(gene_set)
  total_genes <- length(HumanNet_genes) # If you want to enrich against a background of known cancer driver genes instead of against the network
  overlap_count <- length(intersect(gene_set, drivers_in_network))  # Number of overlapping genes
  
  stat <- list()
  
  # Calculate expected overlap by chance
  expected_overlap <- length(drivers_in_network) * (length(gene_set) / total_genes)
  
  # Calculate fold change
  stat$fold_change <- overlap_count / expected_overlap
  
  # Calculate hypergeometric test p-value
  stat$test <- phyper(overlap_count - 1, length(drivers_in_network), total_genes - length(drivers_in_network), length(gene_set), lower.tail = FALSE)
  
  print(paste0(out, ': ', stat$test))
  
  pdf(paste0('../results/scaled/limma/', out, '_brain_Volcano_Abs_with_P.pdf'))
  print(EnhancedVolcano(res,
                        lab = rownames(res),
                        x = 'logFC',
                        y = 'adj.P.Val',
                        labSize = 3,
                        title = out,
                        subtitle = stat$test,
                        pCutoff = 0.05,
                        FCcutoff = 1))
  dev.off()
  #return(stat)
}

test_enrichment(Lioness_lung_diff_nodes, lung_drivers, total_lung, out='LIONESS')
test_enrichment(SSN_lung_diff_nodes, lung_drivers, total_lung, out='SSN')
test_enrichment(CSN_lung_diff_nodes, lung_drivers, total_lung, out='CSN')
test_enrichment(iENA_lung_diff_nodes, lung_drivers, total_lung, out='iENA')
test_enrichment(SSPGI_lung_diff_nodes, lung_drivers, total_lung, out='SSPGI')

test_enrichment(Lioness_brain_diff_nodes, brain_drivers, total_brain, out='LIONESS', lung=FALSE)
test_enrichment(SSN_brain_diff_nodes, brain_drivers, total_brain, out='SSN', lung=FALSE)
test_enrichment(CSN_brain_diff_nodes, brain_drivers, total_brain, out='CSN', lung=FALSE)
test_enrichment(iENA_brain_diff_nodes, brain_drivers, total_brain, out='iENA', lung=FALSE)
test_enrichment(SSPGI_brain_diff_nodes, brain_drivers, total_brain, out='SSPGI', lung=FALSE)


  

