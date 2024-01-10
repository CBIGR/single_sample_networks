#!/usr/bin/Rscript


################################################################################################################
#### Top n connected genes as hubs + check enrichment for driver genes
################################################################################################################

library(plyr)
library(igraph)
library(data.table)
library(UpSetR)
library(reshape2)
library(ggplot2)
library(stringr)
library(ggvenn)
library(ggpubr)
library(stringr)

setwd('/home/boris/Documents/PhD/single_sample/networks')
base <- '/home/boris/Documents/PhD/single_sample/networks/'


convertToSymbol <- function(top200list){
  top200list <- str_replace(top200list, "\\s.*","")
  top200list <- str_replace(top200list, "\\(.*", "")
  return(top200list)
}

top200nodes <- function(network, number_of_drivers, offset=TRUE){
  # network <- paste0(base, "SSPGI/SSPGI_brain_edges_2_75_HumanNet_top_25000_edges.csv")
  # offset <- TRUE
  # number_of_drivers <- 5
  network <- fread(network, data.table=FALSE, header=TRUE, fill=TRUE)
  
  if (offset == TRUE){
    network[,1] <- NULL
  }
  n <- dim(network)[2]-2
  net <- network
  hub_list <- vector("list", length(colnames(network))-2)
  names(hub_list) <- colnames(network)[c(-1,-2)]
  
  for (i in 1:n){
    #i <- 2
    non_null_ind <- which(net[i+2] !=0)
    sample_netw <- net[non_null_ind, c(1,2,i+2)]
    connections <- as.data.frame(sort(table(c(sample_netw$reg, sample_netw$tar)), decreasing = TRUE))
    limit <- connections[number_of_drivers,]$Freq
    top200 <- connections[connections$Freq >= limit, ]
    top200 <- as.vector(top200$Var1)
    top200 <- convertToSymbol(top200)
    hub_list[[colnames(network)[i+2]]] <- top200
  }
  return(hub_list)
}

k <- 200

LIONESS_brain <- top200nodes(paste0(base, "LIONESS/JDS_LIONESS_PCC_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), k, offset=TRUE)
CSN_brain <- top200nodes(paste0(base, "CSN/JDS_CSN_2_75_brain_no_doubles_no_loops_HumanNet.csv"), k, offset=FALSE)
SSN_brain <- top200nodes(paste0(base,"SSN/JDS_SSN_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), k, offset=TRUE)
iENA_brain <- top200nodes(paste0(base,"ssPCC/JDS_iENA_2_75_brain_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), k, offset=TRUE)
SSPGI_brain <- top200nodes(paste0(base,"SSPGI/SSPGI_brain_edges_2_75_HumanNet_top_25000_edges_symbol.csv"), k, offset=TRUE)
SWEET_brain <- top200nodes(paste0(base,"SWEET/SWEET_brain_no_doubles_no_loops_HN_-1_1scaled_top_25000_edges.csv"), k, offset=FALSE)

names(SSN_brain) <- str_remove(names(SSN_brain), "X")
names(iENA_brain) <- str_remove(names(iENA_brain), "X")
sample_info <- fread(paste0(base, "20Q4_v2_sample_info.csv"), header=TRUE, fill=TRUE, data.table=FALSE)
sample_info$DepMap_ID <- str_remove(sample_info$DepMap_ID, "ACH-00")
sample_info_brain <- sample_info[sample_info$DepMap_ID %in% names(LIONESS_brain),]
identical(sample_info_brain$DepMap_ID, names(LIONESS_brain))
#medullo_sample <- sample_info_brain$DepMap_ID[sample_info_brain$Subtype=="Medulloblastoma"]
#GBM_sample <- sample_info_brain$DepMap_ID[sample_info_brain$Subtype=="Glioblastoma"]


# #  network hubs 
aggregate <- fread(paste0(base, "aggregate/bulk_network_brain_PCC_no_doubles_no_loops_HumanNet_top_25000_edges.csv"), data.table=FALSE)
aggregate[,1] <- NULL


connections_agg <- as.data.frame(sort(table(c(aggregate$reg, aggregate$tar)), decreasing = TRUE))
limit_agg <- connections_agg[200,]$Freq
top200_agg <- connections_agg[connections_agg$Freq >= limit_agg, ]
top200_agg <- as.vector(top200_agg$Var1)
top200_agg <- convertToSymbol(top200_agg)


samples <- sample_info_brain$DepMap_ID



agg.overlap <- function(x) {
  return(sum(x %in% top200_agg))
}
# 
# LIONESS_brain_agg <- unlist(lapply(LIONESS_brain, agg.overlap))
# SSN_brain_agg <- unlist(lapply(SSN_brain, agg.overlap))
# iENA_brain_agg <- unlist(lapply(iENA_brain, agg.overlap))
# CSN_brain_agg <- unlist(lapply(CSN_brain, agg.overlap))
# SSPGI_brain_agg <- unlist(lapply(SSPGI_brain, agg.overlap))
# 
# LIONESS_brain_agg_df <- cbind(LIONESS_brain_agg, rep("LIONESS", length(samples)))
# SSN_brain_agg_df <- cbind(SSN_brain_agg, rep("SSN", length(samples)))
# iENA_brain_agg_df <- cbind(iENA_brain_agg, rep("iENA", length(samples)))
# CSN_brain_agg_df <- cbind(CSN_brain_agg, rep("CSN", length(samples)))
# SSPGI_brain_agg_df <- cbind(SSPGI_brain_agg, rep("SSPGI", length(samples)))
# all_agg_df <- data.frame(rbind(LIONESS_brain_agg_df, SSN_brain_agg_df, iENA_brain_agg_df, CSN_brain_agg_df, SSPGI_brain_agg_df))
# colnames(all_agg_df) <- c("amount", "method")
# all_agg_df$amount <- as.numeric((all_agg_df$amount))
# 
# p<-ggplot(all_agg_df, aes(x=method, y=amount, color=method)) +
#   geom_boxplot()
# p
# 
# LIONESS_all <- tail(c(table(factor(table(unlist(LIONESS_brain[samples]))[names(table(unlist(LIONESS_brain[samples])))%in% top200_agg], levels = 1:length(samples)))),1)
# SSN_all <- tail(c(table(factor(table(unlist(SSN_brain[samples]))[names(table(unlist(SSN_brain[samples])))%in% top200_agg], levels = 1:length(samples)))),1)
# iENA_all <- tail(c(table(factor(table(unlist(iENA_brain[samples]))[names(table(unlist(iENA_brain[samples])))%in% top200_agg], levels = 1:length(samples)))),1)
# CSN_all <- tail(c(table(factor(table(unlist(CSN_brain[samples]))[names(table(unlist(CSN_brain[samples])))%in% top200_agg], levels = 1:length(samples)))),1)
# SSPGI_all <- tail(c(table(factor(table(unlist(SSPGI_brain[samples]))[names(table(unlist(SSPGI_brain[samples])))%in% top200_agg], levels = 1:length(samples)))),1)
# 
# p2 <- ggplot() +
#   # box plot of mtcars (mpg vs cyl)
#   geom_boxplot(data = all_agg_df, 
#                aes(x = method, y= amount), fill='#A4A4A4') +
#   # points of data.frame literal
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# pdf(paste0("C:/Users/joked/OneDrive - UGent/PhD Joke shared/paper sample specific benchmarking/Results/Extra_explorations_paper/hub_overlap/overlap_brain_aggregate_all.pdf"), width=8, height=5)
# plot(p2)
# dev.off()
 

# subtype hubs 
subtypes <- c("Glioblastoma", "Medulloblastoma", "Astrocytoma", "Glioma")
samples_1 <- sample_info_brain$DepMap_ID[sample_info_brain$Subtype==subtypes[1]]
samples_not_1 <- sample_info_brain$DepMap_ID[sample_info_brain$Subtype!=subtypes[1]]
samples_2 <- sample_info_brain$DepMap_ID[sample_info_brain$Subtype==subtypes[2]]
samples_not_2 <- sample_info_brain$DepMap_ID[sample_info_brain$Subtype!=subtypes[2]]
samples_3 <- sample_info_brain$DepMap_ID[sample_info_brain$Subtype==subtypes[3]]
samples_not_3 <- sample_info_brain$DepMap_ID[sample_info_brain$Subtype!=subtypes[3]]
samples_4 <- sample_info_brain$DepMap_ID[sample_info_brain$Subtype==subtypes[4]]
samples_not_4 <- sample_info_brain$DepMap_ID[sample_info_brain$Subtype!=subtypes[4]]


top_percentage <- 0.75 
samples <- samples_1
samples_not <- samples_not_1

x <- c(1:length(samples))
LIONESS_y <- table(unlist(LIONESS_brain[samples]))
LIONESS_y_df <- cbind(LIONESS_y, rep("LIONESS", length(LIONESS_y)))
SSN_y <- table(unlist(SSN_brain[samples]))
SSN_y_df <- cbind(SSN_y, rep("SSN", length(SSN_y)))
iENA_y <- table(unlist(iENA_brain[samples]))
iENA_y_df <- cbind(iENA_y, rep("iENA", length(iENA_y)))
CSN_y <- table(unlist(CSN_brain[samples]))
CSN_y_df <- cbind(CSN_y, rep("CSN", length(CSN_y)))
SSPGI_y <- table(unlist(SSPGI_brain[samples]))
SSPGI_y_df <- cbind(SSPGI_y, rep("SSPGI", length(SSPGI_y)))
SWEET_y <- table(unlist(SWEET_brain[samples]))
SWEET_y_df <- cbind(SWEET_y, rep("SWEET", length(SWEET_y)))

all_y_df <- data.frame(rbind(LIONESS_y_df, SSN_y_df, iENA_y_df, CSN_y_df, SSPGI_y_df))
colnames(all_y_df) <- c("amount", "method")
all_y_df$amount <- as.numeric((all_y_df$amount))


LIONESS_hubs_1 <- names(table(unlist(LIONESS_brain[samples]))[table(unlist(LIONESS_brain[samples]))>(length(samples)*top_percentage)])
SSN_hubs_1 <- names(table(unlist(SSN_brain[samples]))[table(unlist(SSN_brain[samples]))>(length(samples)*top_percentage)])
iENA_hubs_1 <- names(table(unlist(iENA_brain[samples]))[table(unlist(iENA_brain[samples]))>(length(samples)*top_percentage)])
CSN_hubs_1 <- names(table(unlist(CSN_brain[samples]))[table(unlist(CSN_brain[samples]))>(length(samples)*top_percentage)])
SSPGI_hubs_1 <- names(table(unlist(SSPGI_brain[samples]))[table(unlist(SSPGI_brain[samples]))>(length(samples)*top_percentage)])
SWEET_hubs_1 <- names(table(unlist(SWEET_brain[samples]))[table(unlist(SWEET_brain[samples]))>(length(samples)*top_percentage)])

LIONESS_hubs_2 <- names(table(unlist(LIONESS_brain[samples_not]))[table(unlist(LIONESS_brain[samples_not]))>(length(samples_not)*top_percentage)])
SSN_hubs_2 <- names(table(unlist(SSN_brain[samples_not]))[table(unlist(SSN_brain[samples_not]))>(length(samples_not)*top_percentage)])
iENA_hubs_2 <- names(table(unlist(iENA_brain[samples_not]))[table(unlist(iENA_brain[samples_not]))>(length(samples_not)*top_percentage)])
CSN_hubs_2 <- names(table(unlist(CSN_brain[samples_not]))[table(unlist(CSN_brain[samples_not]))>(length(samples_not)*top_percentage)])
SSPGI_hubs_2 <- names(table(unlist(SSPGI_brain[samples_not]))[table(unlist(SSPGI_brain[samples_not]))>(length(samples_not)*top_percentage)])
SWEET_hubs_2 <- names(table(unlist(SWEET_brain[samples_not]))[table(unlist(SWEET_brain[samples_not]))>(length(samples_not)*top_percentage)])

LIONESS_un_1 <- LIONESS_hubs_1[!LIONESS_hubs_1 %in% LIONESS_hubs_2]
SSN_un_1 <- SSN_hubs_1[!SSN_hubs_1 %in% SSN_hubs_2] 
iENA_un_1 <- iENA_hubs_1[!iENA_hubs_1 %in% iENA_hubs_2] 
CSN_un_1 <- CSN_hubs_1[!CSN_hubs_1 %in% CSN_hubs_2] 
SSPGI_un_1 <- SSPGI_hubs_1[!SSPGI_hubs_1 %in% SSPGI_hubs_2] 
SWEET_un_1 <- SWEET_hubs_1[!SWEET_hubs_1 %in% SWEET_hubs_2]

LIONESS_rec_df <- cbind(LIONESS_y_df, as.numeric(rownames(LIONESS_y_df) %in% LIONESS_un_1))
SSN_rec_df <- cbind(SSN_y_df, as.numeric(rownames(SSN_y_df) %in% SSN_un_1))
iENA_rec_df <- cbind(iENA_y_df, as.numeric(rownames(iENA_y_df) %in% iENA_un_1))
CSN_rec_df <- cbind(CSN_y_df, as.numeric(rownames(CSN_y_df) %in% CSN_un_1))
SSPGI_rec_df <- cbind(SSPGI_y_df, as.numeric(rownames(SSPGI_y_df) %in% SSPGI_un_1))
SWEET_rec_df <- cbind(SWEET_y_df, as.numeric(rownames(SWEET_y_df) %in% SWEET_un_1))

all_rec_df <- data.frame(rbind(LIONESS_rec_df, SSN_rec_df, SWEET_rec_df, iENA_rec_df, CSN_rec_df, SSPGI_rec_df))
colnames(all_rec_df) <- c("amount", "method", "recurrent")
all_rec_df$amount <- as.numeric((all_rec_df$amount))

result <- all_rec_df %>%
  group_by(method, recurrent) %>%
  summarize(average_amount = mean(amount))

#pdf(paste0("C:/Users/joked/OneDrive - UGent/PhD Joke shared/paper sample specific benchmarking/Results/Extra_explorations_paper/hub_overlap/overlap_brain_", i, ".pdf"), width=8, height=5)
p <- ggplot(all_rec_df, aes(x=method, y=amount)) + 
  geom_violin(trim=FALSE) + geom_jitter( position=position_jitter(0.2), aes(color=recurrent, shape=recurrent)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16), axis.title=element_text(size=16)) +
  xlab("") + ylab("Number of hubs") + 
  scale_color_manual(breaks = c("0", "1"),name="",
                     values=c("black", "#E69F00"),
                     labels = c("Non subtype-specific", "Subtype-specific (75%)"))  +
  scale_shape_manual(breaks = c("0", "1"),name="",
                     values=c(1,0),
                     labels = c("Non subtype-specific", "Subtype-specific (75%)"),
  )
#dev.off()
pdf(paste0("../results/Rebuttal/SWEET/top200_connected_nodes/Glioblastoma_hub_specificity_withSWEET.pdf"), width=8, height=5)
plot(p)
dev.off()

sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="LIONESS"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="SSN"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="iENA"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="CSN"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="SSPGI"]))

#MEDULLOBLASTOMA 

samples <- samples_2
samples_not <- samples_not_2

x <- c(1:length(samples))
LIONESS_y <- table(unlist(LIONESS_brain[samples]))
LIONESS_y_df <- cbind(LIONESS_y, rep("LIONESS", length(LIONESS_y)))
SSN_y <- table(unlist(SSN_brain[samples]))
SSN_y_df <- cbind(SSN_y, rep("SSN", length(SSN_y)))
iENA_y <- table(unlist(iENA_brain[samples]))
iENA_y_df <- cbind(iENA_y, rep("iENA", length(iENA_y)))
CSN_y <- table(unlist(CSN_brain[samples]))
CSN_y_df <- cbind(CSN_y, rep("CSN", length(CSN_y)))
SSPGI_y <- table(unlist(SSPGI_brain[samples]))
SSPGI_y_df <- cbind(SSPGI_y, rep("SSPGI", length(SSPGI_y)))
SWEET_y <- table(unlist(SWEET_brain[samples]))
SWEET_y_df <- cbind(SWEET_y, rep("SWEET", length(SWEET_y)))

all_y_df <- data.frame(rbind(LIONESS_y_df, SSN_y_df, iENA_y_df, CSN_y_df, SSPGI_y_df))
colnames(all_y_df) <- c("amount", "method")
all_y_df$amount <- as.numeric((all_y_df$amount))


LIONESS_hubs_1 <- names(table(unlist(LIONESS_brain[samples]))[table(unlist(LIONESS_brain[samples]))>(length(samples)*top_percentage)])
SSN_hubs_1 <- names(table(unlist(SSN_brain[samples]))[table(unlist(SSN_brain[samples]))>(length(samples)*top_percentage)])
iENA_hubs_1 <- names(table(unlist(iENA_brain[samples]))[table(unlist(iENA_brain[samples]))>(length(samples)*top_percentage)])
CSN_hubs_1 <- names(table(unlist(CSN_brain[samples]))[table(unlist(CSN_brain[samples]))>(length(samples)*top_percentage)])
SSPGI_hubs_1 <- names(table(unlist(SSPGI_brain[samples]))[table(unlist(SSPGI_brain[samples]))>(length(samples)*top_percentage)])
SWEET_hubs_1 <- names(table(unlist(SWEET_brain[samples]))[table(unlist(SWEET_brain[samples]))>(length(samples)*top_percentage)])

LIONESS_hubs_2 <- names(table(unlist(LIONESS_brain[samples_not]))[table(unlist(LIONESS_brain[samples_not]))>(length(samples_not)*top_percentage)])
SSN_hubs_2 <- names(table(unlist(SSN_brain[samples_not]))[table(unlist(SSN_brain[samples_not]))>(length(samples_not)*top_percentage)])
iENA_hubs_2 <- names(table(unlist(iENA_brain[samples_not]))[table(unlist(iENA_brain[samples_not]))>(length(samples_not)*top_percentage)])
CSN_hubs_2 <- names(table(unlist(CSN_brain[samples_not]))[table(unlist(CSN_brain[samples_not]))>(length(samples_not)*top_percentage)])
SSPGI_hubs_2 <- names(table(unlist(SSPGI_brain[samples_not]))[table(unlist(SSPGI_brain[samples_not]))>(length(samples_not)*top_percentage)])
SWEET_hubs_2 <- names(table(unlist(SWEET_brain[samples_not]))[table(unlist(SWEET_brain[samples_not]))>(length(samples_not)*top_percentage)])

LIONESS_un_1 <- LIONESS_hubs_1[!LIONESS_hubs_1 %in% LIONESS_hubs_2]
SSN_un_1 <- SSN_hubs_1[!SSN_hubs_1 %in% SSN_hubs_2] 
iENA_un_1 <- iENA_hubs_1[!iENA_hubs_1 %in% iENA_hubs_2] 
CSN_un_1 <- CSN_hubs_1[!CSN_hubs_1 %in% CSN_hubs_2] 
SSPGI_un_1 <- SSPGI_hubs_1[!SSPGI_hubs_1 %in% SSPGI_hubs_2] 
SWEET_un_1 <- SWEET_hubs_1[!SWEET_hubs_1 %in% SWEET_hubs_2]

LIONESS_rec_df <- cbind(LIONESS_y_df, as.numeric(rownames(LIONESS_y_df) %in% LIONESS_un_1))
SSN_rec_df <- cbind(SSN_y_df, as.numeric(rownames(SSN_y_df) %in% SSN_un_1))
iENA_rec_df <- cbind(iENA_y_df, as.numeric(rownames(iENA_y_df) %in% iENA_un_1))
CSN_rec_df <- cbind(CSN_y_df, as.numeric(rownames(CSN_y_df) %in% CSN_un_1))
SSPGI_rec_df <- cbind(SSPGI_y_df, as.numeric(rownames(SSPGI_y_df) %in% SSPGI_un_1))
SWEET_rec_df <- cbind(SWEET_y_df, as.numeric(rownames(SWEET_y_df) %in% SWEET_un_1))

all_rec_df <- data.frame(rbind(LIONESS_rec_df, SSN_rec_df, SWEET_rec_df, iENA_rec_df, CSN_rec_df, SSPGI_rec_df))
colnames(all_rec_df) <- c("amount", "method", "recurrent")
all_rec_df$amount <- as.numeric((all_rec_df$amount))

result <- all_rec_df %>%
  group_by(method, recurrent) %>%
  summarize(average_amount = mean(amount))

#pdf(paste0("C:/Users/joked/OneDrive - UGent/PhD Joke shared/paper sample specific benchmarking/Results/Extra_explorations_paper/hub_overlap/overlap_brain_", i, ".pdf"), width=8, height=5)
p <- ggplot(all_rec_df, aes(x=method, y=amount)) + 
  geom_violin(trim=FALSE) + geom_jitter( position=position_jitter(0.2), aes(color=recurrent, shape=recurrent)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16), axis.title=element_text(size=16)) + 
  xlab("") + ylab("Number of hubs") + 
  scale_color_manual(breaks = c("0", "1"),name="",
                     values=c("black", "#E69F00"),
                     labels = c("Non subtype-specific", "Subtype-specific (75%)"))  +
  scale_shape_manual(breaks = c("0", "1"),name="",
                     values=c(1,0),
                     labels = c("Non subtype-specific", "Subtype-specific (75%)"),
  )
#dev.off()
pdf(paste0("../results/Rebuttal/SWEET/top200_connected_nodes/Medulloblastoma_hub_specificity.pdf"), width=8, height=5)
plot(p)
dev.off()

sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="LIONESS"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="SSN"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="iENA"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="CSN"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="SSPGI"]))


#ASTROCYTOMA 

samples <- samples_3
samples_not <- samples_not_3

x <- c(1:length(samples))
LIONESS_y <- table(unlist(LIONESS_brain[samples]))
LIONESS_y_df <- cbind(LIONESS_y, rep("LIONESS", length(LIONESS_y)))
SSN_y <- table(unlist(SSN_brain[samples]))
SSN_y_df <- cbind(SSN_y, rep("SSN", length(SSN_y)))
iENA_y <- table(unlist(iENA_brain[samples]))
iENA_y_df <- cbind(iENA_y, rep("iENA", length(iENA_y)))
CSN_y <- table(unlist(CSN_brain[samples]))
CSN_y_df <- cbind(CSN_y, rep("CSN", length(CSN_y)))
SSPGI_y <- table(unlist(SSPGI_brain[samples]))
SSPGI_y_df <- cbind(SSPGI_y, rep("SSPGI", length(SSPGI_y)))
SWEET_y <- table(unlist(SWEET_brain[samples]))
SWEET_y_df <- cbind(SWEET_y, rep("SWEET", length(SWEET_y)))

all_y_df <- data.frame(rbind(LIONESS_y_df, SSN_y_df, iENA_y_df, CSN_y_df, SSPGI_y_df))
colnames(all_y_df) <- c("amount", "method")
all_y_df$amount <- as.numeric((all_y_df$amount))


LIONESS_hubs_1 <- names(table(unlist(LIONESS_brain[samples]))[table(unlist(LIONESS_brain[samples]))>(length(samples)*top_percentage)])
SSN_hubs_1 <- names(table(unlist(SSN_brain[samples]))[table(unlist(SSN_brain[samples]))>(length(samples)*top_percentage)])
iENA_hubs_1 <- names(table(unlist(iENA_brain[samples]))[table(unlist(iENA_brain[samples]))>(length(samples)*top_percentage)])
CSN_hubs_1 <- names(table(unlist(CSN_brain[samples]))[table(unlist(CSN_brain[samples]))>(length(samples)*top_percentage)])
SSPGI_hubs_1 <- names(table(unlist(SSPGI_brain[samples]))[table(unlist(SSPGI_brain[samples]))>(length(samples)*top_percentage)])
SWEET_hubs_1 <- names(table(unlist(SWEET_brain[samples]))[table(unlist(SWEET_brain[samples]))>(length(samples)*top_percentage)])

LIONESS_hubs_2 <- names(table(unlist(LIONESS_brain[samples_not]))[table(unlist(LIONESS_brain[samples_not]))>(length(samples_not)*top_percentage)])
SSN_hubs_2 <- names(table(unlist(SSN_brain[samples_not]))[table(unlist(SSN_brain[samples_not]))>(length(samples_not)*top_percentage)])
iENA_hubs_2 <- names(table(unlist(iENA_brain[samples_not]))[table(unlist(iENA_brain[samples_not]))>(length(samples_not)*top_percentage)])
CSN_hubs_2 <- names(table(unlist(CSN_brain[samples_not]))[table(unlist(CSN_brain[samples_not]))>(length(samples_not)*top_percentage)])
SSPGI_hubs_2 <- names(table(unlist(SSPGI_brain[samples_not]))[table(unlist(SSPGI_brain[samples_not]))>(length(samples_not)*top_percentage)])
SWEET_hubs_2 <- names(table(unlist(SWEET_brain[samples_not]))[table(unlist(SWEET_brain[samples_not]))>(length(samples_not)*top_percentage)])

LIONESS_un_1 <- LIONESS_hubs_1[!LIONESS_hubs_1 %in% LIONESS_hubs_2]
SSN_un_1 <- SSN_hubs_1[!SSN_hubs_1 %in% SSN_hubs_2] 
iENA_un_1 <- iENA_hubs_1[!iENA_hubs_1 %in% iENA_hubs_2] 
CSN_un_1 <- CSN_hubs_1[!CSN_hubs_1 %in% CSN_hubs_2] 
SSPGI_un_1 <- SSPGI_hubs_1[!SSPGI_hubs_1 %in% SSPGI_hubs_2] 
SWEET_un_1 <- SWEET_hubs_1[!SWEET_hubs_1 %in% SWEET_hubs_2]

LIONESS_rec_df <- cbind(LIONESS_y_df, as.numeric(rownames(LIONESS_y_df) %in% LIONESS_un_1))
SSN_rec_df <- cbind(SSN_y_df, as.numeric(rownames(SSN_y_df) %in% SSN_un_1))
iENA_rec_df <- cbind(iENA_y_df, as.numeric(rownames(iENA_y_df) %in% iENA_un_1))
CSN_rec_df <- cbind(CSN_y_df, as.numeric(rownames(CSN_y_df) %in% CSN_un_1))
SSPGI_rec_df <- cbind(SSPGI_y_df, as.numeric(rownames(SSPGI_y_df) %in% SSPGI_un_1))
SWEET_rec_df <- cbind(SWEET_y_df, as.numeric(rownames(SWEET_y_df) %in% SWEET_un_1))

all_rec_df <- data.frame(rbind(LIONESS_rec_df, SSN_rec_df, SWEET_rec_df, iENA_rec_df, CSN_rec_df, SSPGI_rec_df))
colnames(all_rec_df) <- c("amount", "method", "recurrent")
all_rec_df$amount <- as.numeric((all_rec_df$amount))

#pdf(paste0("C:/Users/joked/OneDrive - UGent/PhD Joke shared/paper sample specific benchmarking/Results/Extra_explorations_paper/hub_overlap/overlap_brain_", i, ".pdf"), width=8, height=5)
p <- ggplot(all_rec_df, aes(x=method, y=amount)) + 
  geom_violin(trim=FALSE) + geom_jitter( position=position_jitter(0.2), aes(color=recurrent, shape=recurrent)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16), axis.title=element_text(size=16)) + 
  xlab("") + ylab("Number of hubs") + 
  scale_color_manual(breaks = c("0", "1"),name="",
                     values=c("black", "#E69F00"),
                     labels = c("Non subtype-specific", "Subtype-specific (75%)"))  +
  scale_shape_manual(breaks = c("0", "1"),name="",
                     values=c(1,0),
                     labels = c("Non subtype-specific", "Subtype-specific (75%)"),
  )
#dev.off()
pdf(paste0("../results/Rebuttal/SWEET/top200_connected_nodes/Astrocytoma_hub_specificity.pdf"), width=8, height=5)
plot(p)
dev.off()

sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="LIONESS"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="SSN"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="iENA"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="CSN"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="SSPGI"]))


#GLIOMA 

samples <- samples_4
samples_not <- samples_not_4

x <- c(1:length(samples))
LIONESS_y <- table(unlist(LIONESS_brain[samples]))
LIONESS_y_df <- cbind(LIONESS_y, rep("LIONESS", length(LIONESS_y)))
SSN_y <- table(unlist(SSN_brain[samples]))
SSN_y_df <- cbind(SSN_y, rep("SSN", length(SSN_y)))
iENA_y <- table(unlist(iENA_brain[samples]))
iENA_y_df <- cbind(iENA_y, rep("iENA", length(iENA_y)))
CSN_y <- table(unlist(CSN_brain[samples]))
CSN_y_df <- cbind(CSN_y, rep("CSN", length(CSN_y)))
SSPGI_y <- table(unlist(SSPGI_brain[samples]))
SSPGI_y_df <- cbind(SSPGI_y, rep("SSPGI", length(SSPGI_y)))
SWEET_y <- table(unlist(SWEET_brain[samples]))
SWEET_y_df <- cbind(SWEET_y, rep("SWEET", length(SWEET_y)))

all_y_df <- data.frame(rbind(LIONESS_y_df, SSN_y_df, iENA_y_df, CSN_y_df, SSPGI_y_df))
colnames(all_y_df) <- c("amount", "method")
all_y_df$amount <- as.numeric((all_y_df$amount))


LIONESS_hubs_1 <- names(table(unlist(LIONESS_brain[samples]))[table(unlist(LIONESS_brain[samples]))>(length(samples)*top_percentage)])
SSN_hubs_1 <- names(table(unlist(SSN_brain[samples]))[table(unlist(SSN_brain[samples]))>(length(samples)*top_percentage)])
iENA_hubs_1 <- names(table(unlist(iENA_brain[samples]))[table(unlist(iENA_brain[samples]))>(length(samples)*top_percentage)])
CSN_hubs_1 <- names(table(unlist(CSN_brain[samples]))[table(unlist(CSN_brain[samples]))>(length(samples)*top_percentage)])
SSPGI_hubs_1 <- names(table(unlist(SSPGI_brain[samples]))[table(unlist(SSPGI_brain[samples]))>(length(samples)*top_percentage)])
SWEET_hubs_1 <- names(table(unlist(SWEET_brain[samples]))[table(unlist(SWEET_brain[samples]))>(length(samples)*top_percentage)])

LIONESS_hubs_2 <- names(table(unlist(LIONESS_brain[samples_not]))[table(unlist(LIONESS_brain[samples_not]))>(length(samples_not)*top_percentage)])
SSN_hubs_2 <- names(table(unlist(SSN_brain[samples_not]))[table(unlist(SSN_brain[samples_not]))>(length(samples_not)*top_percentage)])
iENA_hubs_2 <- names(table(unlist(iENA_brain[samples_not]))[table(unlist(iENA_brain[samples_not]))>(length(samples_not)*top_percentage)])
CSN_hubs_2 <- names(table(unlist(CSN_brain[samples_not]))[table(unlist(CSN_brain[samples_not]))>(length(samples_not)*top_percentage)])
SSPGI_hubs_2 <- names(table(unlist(SSPGI_brain[samples_not]))[table(unlist(SSPGI_brain[samples_not]))>(length(samples_not)*top_percentage)])
SWEET_hubs_2 <- names(table(unlist(SWEET_brain[samples_not]))[table(unlist(SWEET_brain[samples_not]))>(length(samples_not)*top_percentage)])

LIONESS_un_1 <- LIONESS_hubs_1[!LIONESS_hubs_1 %in% LIONESS_hubs_2]
SSN_un_1 <- SSN_hubs_1[!SSN_hubs_1 %in% SSN_hubs_2] 
iENA_un_1 <- iENA_hubs_1[!iENA_hubs_1 %in% iENA_hubs_2] 
CSN_un_1 <- CSN_hubs_1[!CSN_hubs_1 %in% CSN_hubs_2] 
SSPGI_un_1 <- SSPGI_hubs_1[!SSPGI_hubs_1 %in% SSPGI_hubs_2] 
SWEET_un_1 <- SWEET_hubs_1[!SWEET_hubs_1 %in% SWEET_hubs_2]

LIONESS_rec_df <- cbind(LIONESS_y_df, as.numeric(rownames(LIONESS_y_df) %in% LIONESS_un_1))
SSN_rec_df <- cbind(SSN_y_df, as.numeric(rownames(SSN_y_df) %in% SSN_un_1))
iENA_rec_df <- cbind(iENA_y_df, as.numeric(rownames(iENA_y_df) %in% iENA_un_1))
CSN_rec_df <- cbind(CSN_y_df, as.numeric(rownames(CSN_y_df) %in% CSN_un_1))
SSPGI_rec_df <- cbind(SSPGI_y_df, as.numeric(rownames(SSPGI_y_df) %in% SSPGI_un_1))
SWEET_rec_df <- cbind(SWEET_y_df, as.numeric(rownames(SWEET_y_df) %in% SWEET_un_1))

all_rec_df <- data.frame(rbind(LIONESS_rec_df, SSN_rec_df, SWEET_rec_df, iENA_rec_df, CSN_rec_df, SSPGI_rec_df))
colnames(all_rec_df) <- c("amount", "method", "recurrent")
all_rec_df$amount <- as.numeric((all_rec_df$amount))

#pdf(paste0("C:/Users/joked/OneDrive - UGent/PhD Joke shared/paper sample specific benchmarking/Results/Extra_explorations_paper/hub_overlap/overlap_brain_", i, ".pdf"), width=8, height=5)
p <- ggplot(all_rec_df, aes(x=method, y=amount)) + 
  geom_violin(trim=FALSE) + geom_jitter( position=position_jitter(0.2), aes(color=recurrent, shape=recurrent)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16), axis.title=element_text(size=16)) + 
  xlab("") + ylab("Number of hubs") + 
  scale_color_manual(breaks = c("0", "1"),name="",
                     values=c("black", "#E69F00"),
                     labels = c("Non subtype-specific", "Subtype-specific (75%)"))  +
  scale_shape_manual(breaks = c("0", "1"),name="",
                     values=c(1,0),
                     labels = c("Non subtype-specific", "Subtype-specific (75%)"),
  )
#dev.off()
pdf(paste0("../results/Rebuttal/SWEET/top200_connected_nodes/Glioma_hub_specificity.pdf"), width=8, height=5)
plot(p)
dev.off()

sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="LIONESS"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="SSN"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="iENA"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="CSN"]))
sum(as.numeric(all_rec_df$recurrent[all_rec_df$method=="SSPGI"]))




# ALL HUBS

samples <- sample_info_brain$DepMap_ID

x <- c(1:length(samples))
LIONESS_y <- table(unlist(LIONESS_brain[samples]))
LIONESS_y_df <- cbind(LIONESS_y, rep("LIONESS", length(LIONESS_y)))
SSN_y <- table(unlist(SSN_brain[samples]))
SSN_y_df <- cbind(SSN_y, rep("SSN", length(SSN_y)))
iENA_y <- table(unlist(iENA_brain[samples]))
iENA_y_df <- cbind(iENA_y, rep("iENA", length(iENA_y)))
CSN_y <- table(unlist(CSN_brain[samples]))
CSN_y_df <- cbind(CSN_y, rep("CSN", length(CSN_y)))
SSPGI_y <- table(unlist(SSPGI_brain[samples]))
SSPGI_y_df <- cbind(SSPGI_y, rep("SSPGI", length(SSPGI_y)))
SWEET_y <- table(unlist(SWEET_brain[samples]))
SWEET_y_df <- cbind(SWEET_y, rep("SWEET", length(SWEET_y)))

all_y_df <- data.frame(rbind(LIONESS_y_df, SSN_y_df, SWEET_y_df, iENA_y_df, CSN_y_df, SSPGI_y_df))
colnames(all_y_df) <- c("amount", "method")
all_y_df$amount <- as.numeric((all_y_df$amount))


LIONESS_rec_df <- cbind(LIONESS_y_df, as.numeric(rownames(LIONESS_y_df) %in% top200_agg))
SSN_rec_df <- cbind(SSN_y_df, as.numeric(rownames(SSN_y_df) %in% top200_agg))
iENA_rec_df <- cbind(iENA_y_df, as.numeric(rownames(iENA_y_df) %in% top200_agg))
CSN_rec_df <- cbind(CSN_y_df, as.numeric(rownames(CSN_y_df) %in% top200_agg))
SSPGI_rec_df <- cbind(SSPGI_y_df, as.numeric(rownames(SSPGI_y_df) %in% top200_agg))
SWEET_rec_df <- cbind(SWEET_y_df, as.numeric(rownames(SWEET_y_df) %in% top200_agg))

all_rec_df <- data.frame(rbind(LIONESS_rec_df, SSN_rec_df, SWEET_rec_df, iENA_rec_df, CSN_rec_df, SSPGI_rec_df))

colnames(all_rec_df) <- c("amount", "method", "recurrent")
all_rec_df$amount <- as.numeric((all_rec_df$amount))

#pdf(paste0("C:/Users/joked/OneDrive - UGent/PhD Joke shared/paper sample specific benchmarking/Results/Extra_explorations_paper/hub_overlap/overlap_brain_", i, ".pdf"), width=8, height=5)
p <- ggplot(all_rec_df, aes(x=method, y=amount)) +
  geom_violin(trim=FALSE) + geom_jitter( position=position_jitter(0.2), aes(color=recurrent, shape=recurrent)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16), axis.title=element_text(size=16)) +
  xlab("") + ylab("Number of hubs") +
  scale_color_manual(breaks = c("0", "1"),name="",
                     values=c("black", "dodgerblue2"),
                     labels = c("Non aggregate hubs", "Aggregate hubs"))  +
  scale_shape_manual(breaks = c("0", "1"),name="",
                     values=c(1,0),
                     labels = c("Non aggregate hubs", "Aggregate hubs"),
  )
#dev.off()
pdf(paste0("../results/Rebuttal/SWEET/top200_connected_nodes/brain_recurrence_aggregate_hubs_overlap.pdf"), width=8, height=5)
plot(p)
dev.off()


