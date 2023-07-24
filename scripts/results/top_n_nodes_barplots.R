#!/usr/bin/Rscript

setwd('/home/boris/Documents/PhD/single_sample/networks')

load('../results/top200_connected_nodes/lung/results_top200.RData')
load('../results/top200_connected_nodes/brain/results200.RData')


# We only need the tables with p-values
# Create a vector of cancer types
cancer_types <- c("NSCLC", "SCLC", "HGG", "MBL")

# Create a vector of tool names
tool_names <- rownames(p_values_NSCLC) # Assuming p_values_NSCLC dataframe has the same row names for all dataframes
tool_names <- str_replace(tool_names, '_lung','')

# Create a data frame to store the 1-p_values for all tools and cancer types
barplot_data <- data.frame(
  Tool = tool_names,
  # NSCLC = p_values_NSCLC[, 1],
  # SCLC = p_values_SCLC[, 1],
  # HGG = p_values_HGG[, 1],
  # MBL = p_values_MBL[, 1]
  # 
  NSCLC = fold_changes_NSCLC[, 1],
  SCLC = fold_changes_SCLC[, 1],
  HGG = fold_changes_HGG[, 1],
  MBL = fold_changes_MBL[, 1]
  
)

# Set the order of cancer types for plotting
#barplot_data <- barplot_data[, order(cancer_types)]

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2") # Colorblind palette
pdf('/home/boris/Documents/PhD/single_sample/results/top200_connected_nodes/barplot_FC_enrich_vs_all_HumanNetGenes.pdf')
# Plot the barplot
barplot(
  as.matrix(barplot_data[, -1]), # Transpose the matrix to have cancer types on the x-axis
  beside = TRUE,
  col = cbbPalette,
  xlab = "Cancer Type",
  ylab = "Fold change",
  #legend.text = barplot_data$Tool,
  #args.legend = list(x = "topright")
)
dev.off()


###########################################################################################################
#### Generate the plot for enrichment vs background of all cancer drivers
###########################################################################################################



setwd('/home/boris/Documents/PhD/single_sample/networks')

load('../results/top200_connected_nodes/lung/results_top200_vs_allDrivers.RData')
load('../results/top200_connected_nodes/brain/results200_vs_allDrivers.RData')


# We only need the tables with p-values
# Create a vector of cancer types
cancer_types <- c("NSCLC", "SCLC", "HGG", "MBL")

# Create a vector of tool names
tool_names <- rownames(p_values_NSCLC) # Assuming p_values_NSCLC dataframe has the same row names for all dataframes
tool_names <- str_replace(tool_names, '_lung','')

# Create a data frame to store the 1-p_values for all tools and cancer types
barplot_data <- data.frame(
  Tool = tool_names,
  # NSCLC = 1 - p_values_NSCLC[, 1],
  # SCLC = 1 - p_values_SCLC[, 1],
  # HGG = 1 - p_values_HGG[, 1],
  # MBL = 1 - p_values_MBL[, 1]
  
  NSCLC = fold_changes_NSCLC[, 1],
  SCLC = fold_changes_SCLC[, 1],
  HGG = fold_changes_HGG[, 1],
  MBL = fold_changes_MBL[, 1]
  
)

# 
# barplot_data <- data.frame(
#   Tool = tool_names,
#   # NSCLC = 1 - p_values_NSCLC[, 1],
#   # SCLC = 1 - p_values_SCLC[, 1],
#   # HGG = 1 - p_values_HGG[, 1],
#   # MBL = 1 - p_values_MBL[, 1]
#   
#   NSCLC = p_values_NSCLC[, 1],
#   SCLC = p_values_SCLC[, 1],
#   HGG = p_values_HGG[, 1],
#   MBL = p_values_MBL[, 1]
#   
# )

# Set the order of cancer types for plotting
#barplot_data <- barplot_data[, order(cancer_types)]

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2") # Colorblind palette

pdf('/home/boris/Documents/PhD/single_sample/results/top200_connected_nodes/barplot_FC_enrich_vs_all_drivers.pdf')
# Plot the barplot
barplot(
  as.matrix(barplot_data[, -1]), # Transpose the matrix to have cancer types on the x-axis
  beside = TRUE,
  col = cbbPalette,
  xlab = "Cancer Type",
  ylab = "Fold change",
  legend.text = barplot_data$Tool,
  args.legend = list(x = "topleft")
)
dev.off()

