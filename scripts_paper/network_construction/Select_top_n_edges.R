#!/usr/bin/Rscript


# Problem: due to the binary nature of edge weights in CSN networks, each sample specific network has less edges compared to those in other methods
# So, find out the average number of edges in these CSN networks, and select the top n edges in other networks to make them comparable in size
setwd('/home/boris/Documents/PhD/single_sample/networks')
CSN_lung <- fread('./CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled.csv', header = TRUE, fill = TRUE, data.table = FALSE)
mean(colSums(CSN_lung[, -c(1,2)])) # 27 813.69
CSN_brain <- fread('./CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled.csv', header = TRUE, fill = TRUE, data.table = FALSE)
mean(colSums(CSN_brain[, -c(1,2)])) # 24 399.31

# Select the top25k edges from LIONESS, SSN, iENA and SSPGI networks in both lung and brain

#' 
#'  
#' 
#' @param n the amount of edges to be selected from the networks 
#' @param table table (for all samples with one method), the first two columns of this table need to be the node names 
#' @param remove_zero_rows if FALSE (default) the whole output table is returned , if TRUE the rows with all zeros are removed from the table 
#' 
#' 
#'

library(data.table)
library(plyr)
library(stringr)


select_top_edges <- function(n, table, remove_zero_rows=FALSE, offset=FALSE){
  #print(dim(table)[2])
  #table <- './LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled'
  # offset <- FALSE
  # remove_zero_rows <- TRUE
  # n = 25000
  out <- paste0(table, '_top_', n, '_edges.csv')
  table <- paste0(table, '.csv')
  table <- fread(table, header = TRUE, fill = TRUE, data.table = FALSE)
  
  if (offset==TRUE){
    table <- table[,-1]
  }
  
  for (i in 3:(dim(table)[2])){
    #print(i) #control
    
    #print(table[1:10,i])
    #print(sort(abs(as.numeric(table[,i])), index.return=TRUE, decreasing=TRUE, na.last=TRUE)$ix[1:10])
    #print(sort(abs(as.numeric(table[,i])), index.return=TRUE, decreasing=TRUE, na.last=TRUE)$ix[(n+1):length(table[,i])])
    
    table[sort(abs(as.numeric(table[,i])), index.return=TRUE, decreasing=TRUE, na.last=TRUE)$ix[(n+1):length(table[,i])],i] <- 0
    
    #the edges are sorted and all non selected edges are set equal to zero
  }
  if (remove_zero_rows == FALSE) {
    write.table(table, out, row.names = FALSE, quote = FALSE, sep = '\t')
  } else if (dim(table)[2]>3){
    table <- table[apply(table[,c(-1,-2)], 1, function(x) !all(x==0)),]
    write.table(table, out, row.names = FALSE, quote = FALSE, sep = '\t')
  } else { # dim table == 3
    table <- table[table[,3]!=0,]
    write.table(table, out, row.names = FALSE, quote = FALSE, sep = '\t')
  }
  
}

setwd('/home/boris/Documents/PhD/single_sample/networks')
select_top_edges(25000, './LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled', remove_zero_rows = TRUE, offset = FALSE)
select_top_edges(25000, './LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled', remove_zero_rows = TRUE, offset = FALSE)
select_top_edges(25000, './SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled', remove_zero_rows = TRUE, offset = FALSE)
select_top_edges(25000, './SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled', remove_zero_rows = TRUE, offset = FALSE)
# select_top_edges(25000, './CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled', remove_zero_rows = TRUE, offset = FALSE)
# select_top_edges(25000, './CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled', remove_zero_rows = TRUE, offset = FALSE)
select_top_edges(25000, './ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet_-1_1scaled', remove_zero_rows = TRUE, offset = FALSE)
select_top_edges(25000, './ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet_-1_1scaled', remove_zero_rows = TRUE, offset = FALSE)
select_top_edges(25000, './SSPGI/SSPGI_brain_edges_2_75_HumanNet_-1_1scaled', remove_zero_rows = TRUE, offset = FALSE) # For SSPGI you need to run the function without writing row.names
select_top_edges(25000, './SSPGI/SSPGI_lung_edges_2_75_HumanNet_-1_1scaled', remove_zero_rows = TRUE, offset = FALSE)
select_top_edges(25000, './aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet_-1_1scaled', remove_zero_rows = TRUE, offset = FALSE)
select_top_edges(25000, './aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet_-1_1scaled', remove_zero_rows = TRUE, offset = FALSE)


##############################################################################################################################
#### Change ensembl ID in SSPGI networks to gene symbol
##############################################################################################################################

# Problem: SSPGI networks contain ensembl ID's instead of gene symbols, this is not compatible with downstream scripts
SSPGI_lung <- fread('./SSPGI/SSPGI_lung_edges_2_75_HumanNet_-1_1scaled_top_25000_edges.csv', header=TRUE, fill = TRUE, data.table = FALSE)
aggregate_lung <- fread('./aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv', header=TRUE, fill = TRUE, data.table = FALSE)
aggregate_lung <- data.frame(unique(c(aggregate_lung[,1], aggregate_lung[,2]))); colnames(aggregate_lung) <- c('ensembl')
rownames(aggregate_lung) <- str_replace(aggregate_lung$ensembl, "\\(.*", "")
aggregate_lung$ensembl <- str_replace(aggregate_lung$ensembl, ".*\\(", ""); aggregate_lung$ensembl <- substr(aggregate_lung$ensembl, 1, nchar(aggregate_lung$ensembl)-1)
symbols <- as.vector(row.names(aggregate_lung)); ensembl <- as.vector(aggregate_lung$ensembl)
reg_symbol <- mapvalues(SSPGI_lung$node1, from = ensembl, to = symbols, warn_missing = FALSE); tar_symbol <- mapvalues(SSPGI_lung$node2, from = ensembl, to = symbols, warn_missing = FALSE)
SSPGI_lung$node1 <- reg_symbol; SSPGI_lung$node2 <- tar_symbol
colnames(SSPGI_lung)[1] <- 'reg'; colnames(SSPGI_lung)[2] <- 'tar'
write.table(SSPGI_lung, './SSPGI/SSPGI_lung_edges_2_75_HumanNet_-1_1scaled_top_25000_edges_symbol.csv', quote =FALSE, sep = ',', row.names=FALSE)

SSPGI_brain <- fread('./SSPGI/SSPGI_brain_edges_2_75_HumanNet_-1_1scaled_top_25000_edges.csv', header=TRUE, fill = TRUE, data.table = FALSE)
aggregate_brain <- fread('./aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet.csv', header=TRUE, fill = TRUE, data.table = FALSE)
aggregate_brain <- data.frame(unique(c(aggregate_brain[,1], aggregate_brain[,2]))); colnames(aggregate_brain) <- c('ensembl')
rownames(aggregate_brain) <- str_replace(aggregate_brain$ensembl, "\\(.*", "")
aggregate_brain$ensembl <- str_replace(aggregate_brain$ensembl, ".*\\(", ""); aggregate_brain$ensembl <- substr(aggregate_brain$ensembl, 1, nchar(aggregate_brain$ensembl)-1)
symbols <- as.vector(row.names(aggregate_brain)); ensembl <- as.vector(aggregate_brain$ensembl)
reg_symbol <- mapvalues(SSPGI_brain$node1, from = ensembl, to = symbols, warn_missing = FALSE); tar_symbol <- mapvalues(SSPGI_brain$node2, from = ensembl, to = symbols, warn_missing = FALSE)
SSPGI_brain$node1 <- reg_symbol; SSPGI_brain$node2 <- tar_symbol
colnames(SSPGI_brain)[1] <- 'reg'; colnames(SSPGI_brain)[2] <- 'tar'
write.table(SSPGI_brain, './SSPGI/SSPGI_brain_edges_2_75_HumanNet_-1_1scaled_top_25000_edges_symbol.csv', quote =FALSE, sep = ',', row.names=FALSE)
