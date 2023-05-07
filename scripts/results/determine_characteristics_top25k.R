#!/usr/bin/Rscipt

#### Script to calculate network characteristics of the single sample networks ####

## Start with aggregate networks ##

library(igraph)
library(data.table)
library(reshape)
library(reshape2)
library(ggplot2)
library(dils)

# For local run
base <- '/home/boris/Documents/PhD/single_sample/networks/'
setwd(paste0(base))

# For HPC:
# base <- '/user/gent/435/vsc43501/single_sample/network_outputs/'
# setwd(base)

#load("top50k_characteristics.RData")
determineCharacteristics <- function(single_sample_network, output_file, offset=TRUE){
  
  # leave offset TRUE when working with top50k networks, use FALSE when working with HumanNet
  # Because of the slightly different file format
  # single_sample_network <- paste0(base, "CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv")
  # offset <- FALSE
  # output_file <- 'CSN_brain_HumanNet_top25k.txt'
  network <- fread(single_sample_network, header=TRUE, data.table=FALSE, fill=TRUE)
  #offset <- FALSE
  if (offset==TRUE){
    network <- network[,-c(1)]
  }
  
  #network <- fread(paste0(base, "LIONESS/JDS_LIONESS_PCC_brain_no_doubles_no_loops_2_75_top_50000_edges.csv"), header=TRUE, data.table=FALSE, fill=TRUE)
  #network <- network[,-c(1)]
  
  network_characteristics <- data.frame(matrix(ncol=6, nrow = ncol(network) - 2))
  colnames <- c("ClusteringCoefficient","density", "NodeBetweenness", "EdgeBetweenness", "diameter", "connected_components") 
  colnames(network_characteristics) <- colnames
  rownames(network_characteristics) <- colnames(network[,-c(1,2)])
  
  n <- dim(network)[2]-2
  net <- network
  
  for (i in 1:n){
    #i <- 13   # For testing purposes
    if (offset==TRUE){
      non_null_ind <- which(net[i+2] !=0) # Extract edges with non-zero weights for a given sample
      sample_netw <- net[non_null_ind, c(1,2,i+2)]
    } else {
      non_null_ind <- which(net[i+2] != 0)
      sample_netw <- net[non_null_ind,c(1,2,i+2)]
    }
    
    # Create edge list with weights for a given sample
    colnames(sample_netw) <- c("From", "To", "Weight")
    single_sample_graph <- graph.data.frame(sample_netw, directed = FALSE)
    single_sample_graph$weight <- sample_netw[,3]
    
    # Clustering coefficient
    clustCoef <- transitivity(single_sample_graph, type="average")
    network_characteristics[i,"ClusteringCoefficient"] <- clustCoef
    
    # Edge density = ratio of the number of edges and the number of possible edges.
    density <- edge_density(single_sample_graph, loops=FALSE) # Considering graphs without loops, since loops were removed in the input graphs
    network_characteristics[i, "density"] <- density
    
    # Node betweenness = the number of shortest paths (geodesics) through a node
    betweenness <- estimate_betweenness(single_sample_graph, vids = V(single_sample_graph), directed = FALSE, cutoff = 0)
    # What do we use here? -> Look this up!
    # Using the average value for now
    av_betweenness <- mean(betweenness)
    network_characteristics[i, "NodeBetweenness"] <- av_betweenness
    
    # Edge betweenness <- the number of shortest paths (geodesics) through an edge
    edge_between <- betweenness(single_sample_graph, v = V(single_sample_graph), directed=FALSE, normalized = TRUE)
    av_edge_between <- mean(edge_between)
    network_characteristics[i, "EdgeBetweenness"] <- av_edge_between
    
    # Diameter
    diameter <- diameter(single_sample_graph, directed=FALSE, unconnected=TRUE) # Unconnected=TRUE -> if unconnected, calculate diameter of connected components and return the largest one
    network_characteristics[i, "diameter"] <- diameter
    
    connected <- count_components(single_sample_graph)
    network_characteristics[i, "connected_components"] <- connected
    #comp_distr <- component_distribution(single_sample_graph)
    #plot(comp_distr)
  }
  write.table(network_characteristics, paste0('../results/network_characteristics/', output_file), sep = "\t", quote=FALSE)
  #network_characteristics$sample_ID <- row.names(network_characteristics) # Not used
  return(network_characteristics)
}

### Characteristics of the aggregate networks ####

network_characteristics_lung <- data.frame(matrix(ncol=6, nrow = 1))
colnames <- c("ClusteringCoefficient","density", "NodeBetweenness", "EdgeBetweenness", "diameter", "components")
colnames(network_characteristics_lung) <- colnames
rownames(network_characteristics_lung) <- c("Aggregate_network")

aggregate_lung <- fread(paste0(base, "./aggregate/bulk_network_lung_PCC_no_doubles_no_loops_top_25000_edges.csv"), header=TRUE, data.table=FALSE, fill=TRUE)
aggregate_lung_HumanNet <- fread(paste0(base, "./aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv"), header=TRUE, data.table=FALSE, fill=TRUE)
aggregate_lung[,1] <- NULL # Not necessary when working with HumanNet networks
num_nodes <- length(unique(as.vector(as.matrix(aggregate_lung_HumanNet[,c(1,2)]))))
aggregate_lung$weight <- abs(aggregate_lung$weight)
aggregate_lung_graph <- graph.data.frame(aggregate_lung, directed=FALSE)
network_characteristics_lung[1,"ClusteringCoefficient"] <- transitivity(aggregate_lung_graph, type="average")
network_characteristics_lung[1,"density"] <- edge_density(aggregate_lung_graph, loops=FALSE)
network_characteristics_lung[1,"EdgeBetweenness"] <- mean(betweenness(aggregate_lung_graph, v = V(aggregate_lung_graph), directed=FALSE, normalized = TRUE))
network_characteristics_lung[1,"NodeBetweenness"] <- mean(estimate_betweenness(aggregate_lung_graph, v = V(aggregate_lung_graph), directed=FALSE, cutoff = 0))
network_characteristics_lung[1,"diameter"] <- diameter(aggregate_lung_graph, directed=FALSE, unconnected=TRUE)
network_characteristics_lung[1,"components"] <- count_components(aggregate_lung_graph)

network_characteristics_brain <- data.frame(matrix(ncol=6, nrow = 1))
colnames <- c("ClusteringCoefficient","density", "NodeBetweenness", "EdgeBetweenness", "diameter", "components")
colnames(network_characteristics_brain) <- colnames
rownames(network_characteristics_brain) <- c("Aggregate_network")

aggregate_brain <- fread(paste0(base, "./aggregate/bulk_network_brain_PCC_no_doubles_no_loops_top_25000_edges.csv"), header=TRUE, data.table=FALSE, fill=TRUE)
aggregate_brain_HumanNet <- fread(paste0(base, "./aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet.csv"), header=TRUE, data.table=FALSE, fill=TRUE)
aggregate_brain[,1] <- NULL # Not necessary when working with HumanNet networks
num_nodes <- length(unique(as.vector(as.matrix(aggregate_brain_HumanNet[,c(1,2)]))))
aggregate_brain$weight <- abs(aggregate_brain$weight)
aggregate_brain_graph <- graph.data.frame(aggregate_brain, directed=FALSE)
network_characteristics_brain[1,"ClusteringCoefficient"] <- transitivity(aggregate_brain_graph, type="average")
network_characteristics_brain[1,"density"] <- edge_density(aggregate_brain_graph, loops=FALSE)
network_characteristics_brain[1,"EdgeBetweenness"] <- mean(betweenness(aggregate_brain_graph, v = V(aggregate_brain_graph), directed=FALSE, normalized = TRUE))
network_characteristics_brain[1,"NodeBetweenness"] <- mean(estimate_betweenness(aggregate_brain_graph, v = V(aggregate_brain_graph), directed=FALSE, cutoff = 0))
network_characteristics_brain[1,"diameter"] <- diameter(aggregate_brain_graph, directed=FALSE, unconnected=TRUE)
network_characteristics_brain[1,"components"] <- count_components(aggregate_brain_graph)

#### Now analyze HumanNet networks ####

single_sample_networks <- list(paste0(base, "LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), 
                               paste0(base, "LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv"),
                               paste0(base, "CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv"), 
                               paste0(base, "CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv"), 
                               paste0(base, "SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv"),
                               paste0(base, "SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv"),
                               paste0(base, "ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet_top_25000_edges.csv"),
                               paste0(base, "ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv"),
                               paste0(base, "SSPGI/SSPGI_brain_edges_2_75_HumanNet_top_25000_edges.csv"),
                               paste0(base, "SSPGI/SSPGI_lung_edges_2_75_HumanNet_top_25000_edges.csv")
)

output_files <- list("LIONESS_lung_HumanNet_top25k-test.txt", 
                     "LIONESS_brain_HumanNet_top25k-test.txt",
                     "CSN_brain_HumanNet_top25k-test.txt",
                     "CSN_brain_HumanNet_top25k-test.txt",
                     "SSN_lung_HumanNet_top25k.txt", 
                     "SSN_brain_HumanNet_top25k.txt",
                     "ssPCC_lung_HumanNet_top25k.txt",
                     "ssPCC_brain_HumanNet_top25k.txt",
                     "SSPGI_brain_HumanNet_top25k.txt",
                     "SSPGI_lung_HumanNet_top25k.txt")

offsets <- list(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
library(parallel)
library(MASS)
numcores <- detectCores()

results <- vector("list", 10)
names(results) <- c('lioness_lung', 'lioness_brain', 'CSN_brain', 'CSN_lung', 'SSN_lung', 'SSN_brain', 'iENA_lung', 'iENA_brain', 'SSPGI_brain', 'SSPGI_lung')

results <- mcmapply(FUN = determineCharacteristics,
         single_sample_networks,
         output_files,
         offsets, 
         mc.cores = getOption("mc.cores", numcores-2))

# results <- list('vector', length=10)
# names(results) <- c('lioness_lung', 'lioness_brain', 'CSN_brain', 'CSN_lung', 'SSN_lung', 'SSN_brain', 'iENA_lung', 'iENA_brain', 'SSPGI_brain', 'SSPGI_lung')

save.image("../results/network_characteristics/HumanNet_top25k_characteristics_parallel.RData")

# load("../results/network_characteristics/HumanNet_top25k_characteristics_parallel.RData")

colnames(results) <- c('lioness_lung', 'lioness_brain', 'CSN_brain', 'CSN_lung', 'SSN_lung', 'SSN_brain', 'iENA_lung', 'iENA_brain', 'SSPGI_brain', 'SSPGI_lung')

res_table <- data.frame(matrix(nrow = 6, ncol = 10))
rownames(res_table) <- row.names(results); colnames(res_table) <- c('lioness_lung', 'lioness_brain', 'CSN_brain', 'CSN_lung', 'SSN_lung', 'SSN_brain', 'iENA_lung', 'iENA_brain', 'SSPGI_brain', 'SSPGI_lung')

# Create table with mean values for paper
for (i in 1:10){
  for (k in 1:6){
    # i <- 1
    # k <- 1
    res_table[k, i] <- mean(results[[k, i]])
  }
}

write.table(res_table, '../results/network_characteristics/overview_table.txt', quote=FALSE, sep='\t')
