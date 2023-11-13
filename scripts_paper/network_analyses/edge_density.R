### Script to create edge score density plots from top 50 000 edges single sample networks ###

library(plyr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(dplyr)
library(data.table)
library(stringr)
library(ggpubr)

setwd('/home/boris/Documents/PhD/single_sample/networks')

CreateDensityPlot <- function(input_network, aggregate, output, offset = TRUE, boundary){
  
  aggregate <- paste0(base, "aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv")
  # input_network <- paste0(base, "SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet.csv")
  # output <- "LIONESS_lung_HumanNet"
  # offset <- FALSE
  # boundary <- 3
  title <- str_replace(output, '_.*','')
  aggregate <- fread(aggregate, header=TRUE, data.table=FALSE, fill=TRUE)
  single_sample_raw <- fread(input_network, data.table=FALSE, header=TRUE, fill=TRUE)
  # This single sample file contains all edges that are present in any top50k network -> most edges are zero in each sample!
  
  # Format network input files if necessary (top50k networks contain an additional first column)
  if (offset == TRUE){
    #aggregate[,1] <- NULL
    single_sample <- single_sample_raw[,-c(1,2,3)]
    single_sample[single_sample==0]=NA
  } else {
    single_sample <- single_sample_raw[,-c(1,2)]
    single_sample[single_sample==0]=NA
    #aggregate[,1] <- NULL
  }
  
  single_sample$average <- rowMeans(single_sample, na.rm=TRUE)
  #single_sample$average
  
  aggregate_and_average <- merge(aggregate, single_sample$average, by=0)
  aggregate_and_average <- as.data.frame(aggregate_and_average[,-c(1:3)])
  colnames(aggregate_and_average) <- c("Aggregate", "Mean")
  
  single_sample <- single_sample[,-c(ncol(single_sample))] # Drop column that contains averages
  single_sample[single_sample==0]=NA # Convert weight zero to NA, these edges are not present in this specific single sample network
  
  single_sample_long <- data.frame(All_weights = c(t(single_sample[,-c(ncol(single_sample))])))
  single_sample_long <- na.omit(single_sample_long) # Nothing should be dropped for HumanNet networks
  network_complete <- join(aggregate_and_average, single_sample_long, type="full")
  #lioness_complete$Row.names <- NULL     # Necessary when using merge, but join function is faster
  network.long <- reshape2::melt(network_complete)
  colnames(network.long) <- c('Edge set', 'Weight')
  #lioness.long <- lioness.long[!(lioness.long$value < -10 | lioness.long$value > 10),]
  #network.long <- network.long[!(network.long$value < -6.5 | network.long$value > 6.5),] # Use 6.5 for HumanNet networks (as in Kuijer paper)
  cbbPalette <- c("#D55E00","#000000", "#009E73")
  graph <- ggplot(network.long, aes(x=as.numeric(Weight), color=`Edge set`)) +
    geom_density(aes(linetype=`Edge set`)) + xlab("Edge weight") + ylab("Density") +
    ggtitle(NULL) + xlim(c(-boundary,boundary)) + theme_bw() + theme(plot.title = element_text(hjust = 0.3), axis.title = element_text(size = 6.5)) + 
    scale_colour_manual(values=cbbPalette) + scale_linetype_manual(values=c('solid', 'dashed', 'solid')) + scale_size_manual(values=c(1,5,1))
    
  
  graph
  # 
  # pdf(paste0('/home/boris/Documents/PhD/single_sample/results/density_plots/HumanNet/', output, ".pdf"))
  # plot(graph)
  # dev.off()
  # 
  #rm(graph)
  return(graph)
}

# For HPC
#base <- "/scratch/gent/vo/000/gvo00027/projects/CBIGR/21JDS_single_sample/network_outputs/"

# For local run
base <- '/home/boris/Documents/PhD/single_sample/networks/'
setwd("/home/boris/Documents/PhD/single_sample/networks")
input_network <- paste0(base, "SSPGI/SSPGI_lung_edges_2_75_HumanNet_top_25000_edges.csv")
network <- fread(input_network, data.table=FALSE, header=TRUE, fill=TRUE)

### On HumanNet networks ###

SSN_lung <- CreateDensityPlot(paste0(base, "SSN/JDS_SSN_2_75_lung_no_doubles_no_loops_HumanNet.csv"), paste0(base, "aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv"), "SSN_lung_HumanNet", offset <- FALSE, 0.3)
SSN_brain <- CreateDensityPlot(paste0(base, "SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet.csv"), paste0(base, "aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet.csv"), "SSN_brain_HumanNet", offset <- FALSE, 0.3)


iENA_lung <- CreateDensityPlot(paste0(base, "ssPCC/JDS_iENA_2_75_lung_no_doubles_no_loops_HumanNet.csv"), paste0(base, "aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv"), "iENA_lung_HumanNet", offset <- FALSE, 3)
iENA_brain <- CreateDensityPlot(paste0(base, "ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet.csv"), paste0(base, "aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet.csv"), "iENA_brain_HumanNet", offset <- FALSE, 3)


lioness_lung <- CreateDensityPlot(paste0(base, "LIONESS/JDS_LIONESS_PCC_2_75_lung_no_doubles_no_loops_HumanNet.csv"), paste0(base, "aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv"), "LIONESS_lung_HumanNet", offset <- FALSE, 3)
lioness_brain <- CreateDensityPlot(paste0(base, "LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet.csv"), paste0(base, "aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet.csv"), "LIONESS_brain_HumanNet", offset <- FALSE, 3)


CSN_brain <- CreateDensityPlot(paste0(base, "CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv"), paste0(base, "aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet.csv"), "CSN_brain_HumanNet", offset <- FALSE, 3)
CSN_lung <- CreateDensityPlot(paste0(base, "CSN/JDS_CSN_2_75_lung_no_doubles_no_loops_HumanNet.csv"), paste0(base, "aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv"), "CSN_lung_HumanNet", offset <- FALSE, 3)


SSPGI_brain <- CreateDensityPlot(paste0(base, "SSPGI/SSPGI_brain_edges_2_75_HumanNet.csv"), paste0(base, "aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet.csv"), "SSPGI_lung_HumanNet", offset = TRUE, 5000)
SSPGI_lung <- CreateDensityPlot(paste0(base, "SSPGI/SSPGI_lung_edges_2_75_HumanNet.csv"), paste0(base, "aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv"), "SSPGI_brain_HumanNet", offset=TRUE, 5000)


SWEET_lung <- CreateDensityPlot(paste0(base,"SWEET/SWEET_lung_no_doubles_no_loops_HN.csv"), paste0(base, "aggregate/bulk_network_lung_PCC_no_doubles_no_loops_HumanNet.csv"), "SSPGI_brain_HumanNet", offset=FALSE, 3)
SWEET_brain <- CreateDensityPlot(paste0(base,"SWEET/SWEET_brain_no_doubles_no_loops_HN.csv"), paste0(base, "aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet.csv"), "SSPGI_brain_HumanNet", offset=FALSE, 3)


SWEET_lung_net <- fread(paste0(base,"SWEET/SWEET_lung_no_doubles_no_loops_HN.csv"), header=TRUE, data.table=FALSE, fill=TRUE)
SWEET_brain_net <- fread(paste0(base,"SWEET/SWEET_brain_no_doubles_no_loops_HN.csv"), header=TRUE, data.table=FALSE, fill=TRUE)

# Min weight
min(SWEET_lung_net[, -c(1,2)])
min(SWEET_brain_net[, -c(1,2)])

# Max weight
max(SWEET_lung_net[, -c(1,2)])
max(SWEET_brain_net[, -c(1,2)])


leg <- get_legend(lioness_lung)
lioness_lung <- lioness_lung + theme(legend.position = 'none')
SSN_lung  <- SSN_lung + theme(legend.position = 'none')
CSN_lung <- CSN_lung + theme(legend.position = 'none')
iENA_lung <- iENA_lung + theme(legend.position = 'none')
SSPGI_lung <- SSPGI_lung + theme(legend.position = 'none')
SWEET_lung <- SWEET_lung + theme(legend.position = 'none')

pdf('../results/density_plots/HumanNet/Edge_densities_lung_HumanNet.pdf')
ggarrange(SSN_lung, lioness_lung, iENA_lung, CSN_lung, SSPGI_lung, leg, 
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 3)
dev.off()

leg <- get_legend(lioness_brain)
lioness_brain <- lioness_brain + theme(legend.position = 'none')
SSN_brain  <- SSN_brain + theme(legend.position = 'none')
CSN_brain <- CSN_brain + theme(legend.position = 'none')
iENA_brain <- iENA_brain + theme(legend.position = 'none')
SSPGI_brain <- SSPGI_brain + theme(legend.position = 'none') 
SWEET_brain <- SWEET_brain + theme(legend.position = 'none')

pdf('../results/density_plots/HumanNet/Edge_densities_brain_HumanNet.pdf')
ggarrange(SSN_brain, lioness_brain, iENA_brain, CSN_brain, SSPGI_brain, leg,
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 3)
dev.off()


pdf('../results/Rebuttal/Edge_densities_combination_HumanNet.pdf')
ggarrange(SSN_lung, NULL, SSN_brain, lioness_lung, NULL, lioness_brain,SWEET_lung, NULL, SWEET_brain, iENA_lung, NULL, iENA_brain, CSN_lung, NULL, CSN_brain, SSPGI_lung, NULL, SSPGI_brain,
          # labels = c("A","", "B", "C","", "D", "E","","F","G","", "H", "I","", "J"),
          ncol = 3, nrow = 6, widths = c(1,0.3,1,1,0.3,1,1,0.3,1,1,0.3,1,1,0.3,1,1,0.3,1), common.legend = TRUE)
dev.off()
