#!/usr/bin/Rscript

setwd('/home/boris/Documents/PhD/single_sample/networks')

load('../results/top200_connected_nodes/lung/results_top200.RData')
load('../results/top200_connected_nodes/brain/results200.RData')

# Create a vector of cancer types
cancer_types <- c("NSCLC", "SCLC", "HGG", "MBL")

# Create a vector of tool names
tool_names <- rownames(fold_changes_NSCLC) # Assuming fold_changes_NSCLC dataframe has the same row names for all dataframes
tool_names <- str_replace(tool_names, '_lung', '')

# Create a data frame to store the fold changes and p-values for all tools and cancer types
barplot_data <- data.frame(
  Tool = tool_names,
  NSCLC = fold_changes_NSCLC[, 1],
  SCLC = fold_changes_SCLC[, 1],
  HGG = fold_changes_HGG[, 1],
  MBL = fold_changes_MBL[, 1],
  p_values_NSCLC = p_values_NSCLC[, 1],
  p_values_SCLC = p_values_SCLC[, 1],
  p_values_HGG = p_values_HGG[, 1],
  p_values_MBL = p_values_MBL[, 1]
)

# Set the order of cancer types for plotting
#barplot_data <- barplot_data[, order(cancer_types)]

# Create a color palette for p-values
p_value_colors <- colorRampPalette(c("white", "blue"))(100) # Adjust the color range as desired

# Plot the barplot with colored bars based on p-values
barplot(
  t(as.matrix(barplot_data[, 2:5])), # Transpose the matrix to have cancer types on the x-axis, excluding the p-value columns
  beside = TRUE,
  col = p_value_colors[cut(barplot_data[, 6:9], breaks = 100)],
  xlab = "Cancer Type",
  ylab = "Fold change",
  ylim = c(0, max(barplot_data[, 2:5], na.rm = TRUE)), # Set the y-axis limit to the maximum fold change value
  axis.lty = 1
)

# Add tool name annotations under the bars
text(
  x = barplot_data[, 2:5],
  y = -0.2 * max(barplot_data[, 2:5], na.rm = TRUE), # Adjust the vertical position of the annotations
  labels = barplot_data$Tool,
  xpd = TRUE,
  srt = 45,
  adj = c(1, 0.5),
  cex = 0.8
)
