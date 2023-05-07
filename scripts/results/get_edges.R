#!/usr/bin/Rscipt

library(data.table)
library(stringr)

setwd('/home/boris/Documents/PhD/single_sample/networks/')

LIONESS <- fread('./LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', header=TRUE, data.table=FALSE, fill=TRUE)
SSN <- fread('./SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', header=TRUE, data.table=FALSE, fill=TRUE)
iENA <- fread('./ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', header=TRUE, data.table=FALSE, fill=TRUE)
CSN <- fread('./CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv', header=TRUE, data.table=FALSE, fill=TRUE)
SSPGI <- fread('./SSPGI/SSPGI_brain_edges_2_75_HumanNet_-1_1scaled_top_25000_edges_symbol.csv', header=TRUE, data.table=FALSE, fill=TRUE)  

get_edges <- function(gene, network, folder) {
  # Convert id's to gene symbol
  network <- LIONESS
  gene <- 'EGFR'
  folder = './LIONESS/'
  
  colnames(network) <- gsub('X', '', colnames(network))
  network$reg <- str_replace(network$reg, "\\s.*", ""); network$reg <- str_replace(network$reg, "\\(.*", "")
  network$tar <- str_replace(network$tar, "\\s.*", ""); network$tar <- str_replace(network$tar, "\\(.*", "")
  
  ind <- which(network[,1] %in% gene)
  ind <- c(ind, which(network[,2] %in% gene))
  network <- network[ind, ]
  colSums(network != 0)
  apply(network, 2, function(c)sum(c!=0.000000000)) # -> to many edges per sample to plot -> get the top 10 edges
  
  for (sample in colnames(network[,-c(1,2)])){
    #sample <- "0025"
    to_keep <- as.data.frame(network[, c('reg', 'tar', sample)])
    to_keep[,3] <- abs(to_keep[,3])
    to_keep <- to_keep[order(to_keep[,3], decreasing = TRUE), ]
    to_keep <- to_keep[c(1:10), ]
    to_keep$interaction_type <- 'PP'
    to_keep <- to_keep[, c(1,4,2,3)]
    
    genes <- unique(as.vector(as.matrix(to_keep[,c(1,3)])))
    print(paste(sample, intersect(genes, brain_drivers)))
    
    write.table(to_keep, file = paste0(folder, 'sample_networks/', sample, '_', gene, '.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
    }
}

HGG_drivers <- fread("/home/boris/Documents/PhD/single_sample/networks/IntOGen-DriverGenes_GBM.tsv")[,1]
MBL_drivers <- fread("/home/boris/Documents/PhD/single_sample/networks/IntOGen-DriverGenes_MBL.tsv")[,1]
brain_drivers <- rbind(HGG_drivers, MBL_drivers); brain_drivers <- unique(as.vector(brain_drivers$Symbol))

get_edges('EGFR', LIONESS, './LIONESS/')
get_edges('EGFR', SSN, './SSN/')
get_edges('EGFR', iENA, './ssPCC/')
get_edges('EGFR', CSN, './CSN/')
get_edges('EGFR', SSPGI, './SSPGI/')

#### Get edges connected to EGFR in the aggregate network
aggregate <- fread('./aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet_-1_1scaled_top_25000_edges.csv', header=TRUE, data.table=FALSE, fill=TRUE)
gene <- 'EGFR'; network <- aggregate

colnames(network) <- gsub('X', '', colnames(network))
network$reg <- str_replace(network$reg, "\\s.*", ""); network$reg <- str_replace(network$reg, "\\(.*", "")
network$tar <- str_replace(network$tar, "\\s.*", ""); network$tar <- str_replace(network$tar, "\\(.*", "")

ind <- which(network[,1] %in% gene)
ind <- c(ind, which(network[,2] %in% gene))
to_keep <- network[ind, ]

to_keep[,3] <- abs(to_keep[,3])
to_keep <- to_keep[order(to_keep[,3], decreasing = TRUE), ]
to_keep <- to_keep[c(1:10), ]
to_keep$interaction_type <- 'PP'
to_keep <- to_keep[, c(1,4,2,3)]

genes <- unique(as.vector(as.matrix(to_keep[,c(1,3)])))
print(paste(sample, intersect(genes, brain_drivers)))

write.table(to_keep, file = paste0('./aggregate/', gene, '.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
