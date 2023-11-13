#!/usr/bin/Rscript

library(stringr)

setwd('/home/boris/Documents/PhD/single_sample/networks')



######################################################################################################
#### Main figure: NSCLC, SCLC, MBL and GBM IntOGen drivers
######################################################################################################


load('../results/Rebuttal/top200_connected_nodes/results_top200_lung_Intogen_drivers_withSWEET.RData')
load('../results/Rebuttal/top200_connected_nodes/results_top200_brain_Intogen_drivers_withSWEET.RData')

cancer_types <- c("SCLC", 'NSCLC', 'MBL', 'GBM')
tool_names <- rownames(fold_changes_SCLC) # Assuming fold_changes_NSCLC dataframe has the same row names for all dataframes
tool_names <- str_replace(tool_names, '_lung', '')
tool_names <- c('LIONESS', 'SSN', 'iENA', 'CSN', 'SSPGI', 'SWEET')

# Create a data frame to store the 1-p_values for all tools and cancer types
barplot_data <- data.frame(
  Tool = tool_names,
  NSCLC = fold_changes_NSCLC[, 1],
  SCLC = fold_changes_SCLC[, 1],
  GBM = fold_changes_HGG[, 1],
  MBL = fold_changes_MBL[, 1]
  
)


barplot_data_complete <- data.frame(
  Tool = tool_names,
  NSCLC = fold_changes_NSCLC[, 1],
  SCLC = fold_changes_SCLC[, 1],
  GBM = fold_changes_HGG[, 1],
  MBL = fold_changes_MBL[, 1],
  
  NSCLC_p = p_values_NSCLC[, 1],
  SCLC_p = p_values_SCLC[, 1],
  HGG_p = p_values_HGG[, 1],
  MBL_p = p_values_MBL[, 1]
)
barplot_data <- barplot_data[c(2,1,6,3,4,5), ]
barplot_data_complete <- barplot_data_complete[c(2,1,6,3,4,5), ]
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00") # Colorblind palette
pdf('/home/boris/Documents/PhD/single_sample/results/Rebuttal//SWEET/top200_connected_nodes/barplot_main_IntOGen_drivers.pdf')
# Plot the barplot
barplot(
  as.matrix(barplot_data[, -1]), # Transpose the matrix to have cancer types on the x-axis
  beside = TRUE,
  col = cbbPalette,
  xlab = "Cancer Type",
  ylab = "Fold change",
  ylim = c(0,4),
  #legend.text = barplot_data$Tool,
  #args.legend = list(x = "topleft"),
  main = "IntOGen/COSMIC drivers"
)
dev.off()

######################################################################################################
#### Main figure: NSCLC, SCLC, MBL and GBM Cell Model Passports
######################################################################################################

load('../results/Rebuttal/top200_connected_nodes/results_top200_lung_CellPassports_drivers_withSWEET.RData')
load('../results/Rebuttal/top200_connected_nodes/results_top200_brain_CellPassports_drivers_withSWEET.RData')

cancer_types <- c("SCLC", 'NSCLC', 'MBL', 'GBM')
tool_names <- rownames(fold_changes_SCLC) # Assuming fold_changes_NSCLC dataframe has the same row names for all dataframes
tool_names <- str_replace(tool_names, '_lung', '')
tool_names <- c('LIONESS', 'SSN', 'iENA', 'CSN', 'SSPGI', 'SWEET')

# Create a data frame to store the 1-p_values for all tools and cancer types
barplot_data <- data.frame(
  Tool = tool_names,
  NSCLC = fold_changes_NSCLC[, 1],
  SCLC = fold_changes_SCLC[, 1],
  GBM = fold_changes_HGG[, 1],
  MBL = fold_changes_MBL[, 1]
  
)


barplot_data_complete <- data.frame(
  Tool = tool_names,
  NSCLC = fold_changes_NSCLC[, 1],
  SCLC = fold_changes_SCLC[, 1],
  GBM = fold_changes_HGG[, 1],
  MBL = fold_changes_MBL[, 1],
  
  NSCLC_p = p_values_NSCLC[, 1],
  SCLC_p = p_values_SCLC[, 1],
  HGG_p = p_values_HGG[, 1],
  MBL_p = p_values_MBL[, 1]
)

barplot_data <- barplot_data[c(2,1,6,3,4,5), ]
barplot_data_complete <- barplot_data_complete[c(2,1,6,3,4,5), ]

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00") # Colorblind palette
pdf('/home/boris/Documents/PhD/single_sample/results/Rebuttal/SWEET/top200_connected_nodes/barplot_main_CellPassport_drivers.pdf')
# Plot the barplot
barplot(
  as.matrix(barplot_data[, -1]), # Transpose the matrix to have cancer types on the x-axis
  beside = TRUE,
  col = cbbPalette,
  xlab = "Cancer Type",
  ylab = "Fold change",
  ylim = c(0,15),
  legend.text = barplot_data$Tool,
  args.legend = list(x = "topleft"),
  main = "Cell Model Passports drivers"
)
dev.off()


######################################################################################################
#### Supplementary figure: NSCLC subtypes and SCLC IntOGen drivers
######################################################################################################


load('../results/Rebuttal/top200_connected_nodes/results_top200_lung_Intogen_drivers_withSWEET.RData')

# Create a vector of cancer types
cancer_types <- c("SCLC", "NSCLC_Adenocarcinoma", "NSCLC_Squamous", 'NSCLC')
#cancer_types <- c("SCLC", "NSCLC_Adenocarcinoma", "NSCLC_Squamous")

# Create a vector of tool names
tool_names <- rownames(fold_changes_SCLC) # Assuming fold_changes_NSCLC dataframe has the same row names for all dataframes
tool_names <- str_replace(tool_names, '_lung', '')
tool_names <- c('LIONESS', 'SSN', 'iENA', 'CSN', 'SSPGI', 'SWEET')

# Create a data frame to store the 1-p_values for all tools and cancer types
barplot_data <- data.frame(
  Tool = tool_names,
  # NSCLC = p_values_NSCLC[, 1],
  # SCLC = p_values_SCLC[, 1],
  # HGG = p_values_HGG[, 1],
  # MBL = p_values_MBL[, 1]
  # 
  SCLC = fold_changes_SCLC[, 1],
  NSCLC_Adeno = fold_changes_Adeno[, 1],
  NSCLC_Squamous = fold_changes_Squamous[, 1],
  NSCLC = fold_changes_NSCLC[, 1]
  
)

barplot_data <- barplot_data[c(2,1,6,3,4,5), ]

barplot_data_complete <- data.frame(
  Tool = tool_names,

  SCLC_p = p_values_SCLC[, 1],
  NSCLC_Adeno_p = p_values_Adeno[, 1],
  NSCLC_Squamous_p = p_values_Squamous[, 1],

  NSCLC_p = p_values_NSCLC[, 1],
  #
  SCLC = fold_changes_SCLC[, 1],
  NSCLC_Adeno = fold_changes_Adeno[, 1],
  NSCLC_Squamous = fold_changes_Squamous[, 1],

  NSCLC = fold_changes_NSCLC[, 1]

)

barplot_data_complete <- barplot_data_complete[c(2,1,6,3,4,5), ]

 # Colorblind palette


pdf('/home/boris/Documents/PhD/single_sample/results/Rebuttal/SWEET/top200_connected_nodes/barplot_suppl_lung_IntoGen.pdf')
# Plot the barplot
barplot(
  as.matrix(barplot_data[, -1]), # Transpose the matrix to have cancer types on the x-axis
  beside = TRUE,
  col = cbbPalette,
  xlab = "Cancer Type",
  ylab = "Fold change",
  ylim = c(0,5),
  legend.text = barplot_data$Tool,
  args.legend = list(x = "topleft"),
  main = "IntOGen/COSMIC drivers"
)
dev.off()

######################################################################################################
#### Supplementary figure: NSCLC subtypes and SCLC Cell Model Passports drivers
######################################################################################################

# Same figure both for cellpassport driver genes

load('../results/Rebuttal/top200_connected_nodes/results_top200_lung_CellPassports_drivers_withSWEET.RData')


# We only need the tables with p-values
# Create a vector of cancer types
cancer_types <- c("SCLC", "NSCLC_Adenocarcinoma", "NSCLC_Squamous", "NSCLC_LargeCell", 'NSCLC')
#cancer_types <- c("SCLC", "NSCLC_Adenocarcinoma", "NSCLC_Squamous")

# Create a vector of tool names
tool_names <- rownames(fold_changes_SCLC) # Assuming fold_changes_NSCLC dataframe has the same row names for all dataframes
tool_names <- str_replace(tool_names, '_lung', '')
tool_names <- c('LIONESS', 'SSN', 'iENA', 'CSN', 'SSPGI', 'SWEET')

# Create a data frame to store the 1-p_values for all tools and cancer types
barplot_data <- data.frame(
  Tool = tool_names,
  # NSCLC = p_values_NSCLC[, 1],
  # SCLC = p_values_SCLC[, 1],
  # HGG = p_values_HGG[, 1],
  # MBL = p_values_MBL[, 1]
  # 
  SCLC = fold_changes_SCLC[, 1],
  NSCLC_Adeno = fold_changes_Adeno[, 1],
  NSCLC_Squamous = fold_changes_Squamous[, 1],
  NSCLC_LargeCell = fold_changes_LargeCell[, 1],
  NSCLC = fold_changes_NSCLC[, 1]
)

barplot_data <- barplot_data[c(2,1,6,3,4,5), ]

barplot_data_complete <- data.frame(
  Tool = tool_names,
  
  SCLC_p = p_values_SCLC[, 1],
  NSCLC_Adeno_p = p_values_Adeno[, 1],
  NSCLC_Squamous_p = p_values_Squamous[, 1], 
  NSCLC_LargeCell_p = p_values_LargeCell[, 1],
  NSCLC_p = p_values_NSCLC[, 1],
  # 
  SCLC = fold_changes_SCLC[, 1],
  NSCLC_Adeno = fold_changes_Adeno[, 1],
  NSCLC_Squamous = fold_changes_Squamous[, 1],
  NSCLC_LargeCell = fold_changes_LargeCell[, 1],
  NSCLC = fold_changes_NSCLC[, 1]
  
)
barplot_data_complete <- barplot_data_complete[c(2,1,6,3,4,5), ]

# Set the order of cancer types for plotting
#barplot_data <- barplot_data[, order(cancer_types)]


pdf('/home/boris/Documents/PhD/single_sample/results/Rebuttal/SWEET/top200_connected_nodes/barplot_suppl_lung_CellPassport.pdf')
# Plot the barplot
barplot(
  as.matrix(barplot_data[, -1]), # Transpose the matrix to have cancer types on the x-axis
  beside = TRUE,
  col = cbbPalette,
  xlab = "Cancer Type",
  ylab = "Fold change",
  ylim = c(0,7),
  legend.text = barplot_data$Tool,
  args.legend = list(x = "topleft"),
  main = 'Cell Model Passports drivers'
)
dev.off()

######################################################################################################
#### Supplementary figure: Brain subtypes Cell Model Passports drivers
######################################################################################################

load('../results/Rebuttal/top200_connected_nodes/results_top200_brain_CellPassports_drivers_withSWEET.RData')

cancer_types <- c("GBM", "MBL", "Astrocytoma", 'Oligodendroglioma', 'Glioma')
#cancer_types <- c("SCLC", "NSCLC_Adenocarcinoma", "NSCLC_Squamous")

# Create a vector of tool names
tool_names <- rownames(fold_changes_GBM) # Assuming fold_changes_NSCLC dataframe has the same row names for all dataframes
tool_names <- str_replace(tool_names, '_brain', '')
tool_names <- c('LIONESS', 'SSN', 'iENA', 'CSN', 'SSPGI', 'SWEET')

# Create a data frame to store the 1-p_values for all tools and cancer types
barplot_data <- data.frame(
  Tool = tool_names,
  
  
  GBM = fold_changes_GBM[, 1],
  MBL = fold_changes_MBL[, 1],
  Astrocytoma = fold_changes_astrocytoma[, 1],
  #Oligodendroglioma = fold_changes_oligodendroglioma[, 1],
  #Meningioma = fold_changes_meningioma[, 1],
  Glioma = fold_changes_glioma[, 1]
  
)

barplot_data <- barplot_data[c(2,1,6,3,4,5), ]

barplot_data_complete <- data.frame(
  Tool = tool_names,
  GBM_p = p_values_GBM[, 1],
  MBL_p = p_values_MBL[, 1],
  Astrocytoma_p = p_values_astrocytoma[, 1],
  #Oligodendroglioma = fold_changes_oligodendroglioma[, 1],
  #Meningioma = fold_changes_meningioma[, 1],
  Glioma_p = p_values_glioma[, 1],
  
  GBM = fold_changes_GBM[, 1],
  MBL = fold_changes_MBL[, 1],
  Astrocytoma = fold_changes_astrocytoma[, 1],
  #Oligodendroglioma = fold_changes_oligodendroglioma[, 1],
  #Meningioma = fold_changes_meningioma[, 1],
  Glioma = fold_changes_glioma[, 1]
)

barplot_data_complete <- barplot_data_complete[c(2,1,6,3,4,5), ]

# Set the order of cancer types for plotting
#barplot_data <- barplot_data[, order(cancer_types)]


pdf('/home/boris/Documents/PhD/single_sample/results/Rebuttal/SWEET/top200_connected_nodes/barplot_suppl_brain_CellPassport.pdf')
# Plot the barplot
barplot(
  as.matrix(barplot_data[, -1]), # Transpose the matrix to have cancer types on the x-axis
  beside = TRUE,
  col = cbbPalette,
  xlab = "Cancer Type",
  ylab = "Fold change",
  ylim = c(0,15),
  legend.text = barplot_data$Tool,
  args.legend = list(x = "topleft"),
  main = "Cell Model Passports drivers"
)
dev.off()
