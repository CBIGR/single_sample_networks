#!/usr/bin/Rscript

library("data.table")
library("funkyheatmap")
library(ggplot2)
library(dplyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)

data("dynbenchmark_data")


data <- data.frame(
  "Tools" = c('SSN', 'LIONESS', 'iENA', 'CSN', 'SSPGI', 'SWEET'),
  "Aggregate_network" = c('Yes', 'Yes', 'No', 'No', 'No', 'Yes'),
  "Background_network" = c('No', 'No', 'No', 'No', 'Yes', 'No'),
  "Scalable" = c('Yes', 'Yes', 'Yes', 'Yes', 'No', 'Yes'),
  "Implementation" = c('R', 'R', 'R', 'Matlab', 'R', 'Python'),
  #"Ease-of-use" = c(16, 17, 18, 19, 20),
  
  "Allows_advanced_GRN_inference" = c('No', 'Yes', 'No', 'No', 'No', 'No'),
  "Edge_bias" = c('Yes', 'Yes', 'Yes', 'No', 'No', 'No'),
  
  "Subtype_specificity_NSCLC" = c(57.60000/18.76434, 57.50000/19.37353, 57.60000/19.47798, 61.81818/34.62143, 58.85714/22.88298, 0/35.29291),
  "Subtype_specificity_SCLC" = c(10.350000/3.807947, 10.450000/3.919795,10.478261/3.852349, 11.136364/6.864307, 10.575758/4.652268, 10.800000/8.607914),
  "Subtype_specificity_GBM" = c(28.80000/12.23386, 30.00000/12.90018, 28.88000/12.64935, 29.33333/17.19477, 28.57143/12.42024, 30.00000/21.81579),
  "Subtype_specificity_MBL" = c(7.390244/3.726636, 7.543478/3.397380, 7.407407/3.439359, 7.500000/5.194268, 7.750000/3.744076, 7.363636/7.108000),
  
  "Enrichment_drivers_NSCLC" = 1-c(0.0037694326, 0.0027385331, 0.0028922484, 0.2717082490, 0.0009425231, 0.0079578945),
  "Enrichment_drivers_SCLC" = 1-c(0.1142400, 0.1038216, 0.1112935, 1.0000000, 0.1508792, 0.1278170),
  "Enrichment_drivers_GBM" = 1-c(0.001908186, 0.007868808, 0.001413897, 0.011274148, 0.010208182, 0.024420193),
  "Enrichment_drivers_MBL" = 1-c(0.2469600, 0.2213832, 0.2374178, 0.1360434, 0.2163160, 1.0000000),
  
  "Enrichment_mutations_NSCLC" = 1-c(0.007080154, 0.005145289, 0.005454732, 0.001453609, 0.001405961, 0.001903813),
  "Enrichment_mutations_SCLC" = 1-c(0.131481747, 0.119801139, 0.128182623, 0.070123286, 0.004328369, 0.741135434),
  "Enrichment_mutations_GBM" = 1-c(0.001908186, 0.007868808, 0.001413897, 0.011274148, 0.010208182, 0.024420193),
  "Enrichment_mutations_MBL" = 1-c(0.011281068, 0.009767207, 0.010706050, 0.005311863, 0.009477518, 1.0000000),
  
  "Correlation_lung_proteomics" = c(0.12941561, 0.13751287, 0.07020109, 0.10166485, 0.05580398, 0.10578565),
  "Correlation_brain_proteomics" = c(0.17035069, 0.21115339, 0.17583510, 0.06719328, 0.13434000, 0.21373784),
  "Correlation_lung_CNV" = c(0.03600471, 0.02953851, 0.01968635, 0.01820532, 0.02786514, 0.03320104),
  "Correlation_brain_CNV" = c(0.02858101, 0.02819691, 0.01737868, 0.01831697, 0.02636951, 0.02932895)
  
)




# Rescale the data
rescale1 <- c('Enrichment_drivers_NSCLC', 'Enrichment_drivers_SCLC', 'Enrichment_drivers_GBM', 'Enrichment_drivers_MBL', 'Enrichment_mutations_NSCLC', 'Enrichment_mutations_SCLC', 'Enrichment_mutations_GBM', 'Enrichment_mutations_MBL')
data[, rescale1] <- (data[, rescale1] - 0.5) * 2

rescale2 <- c('Correlation_lung_proteomics', 'Correlation_brain_proteomics', 'Correlation_lung_CNV', 'Correlation_brain_CNV')
data[, rescale2] <- (data[, rescale2]) * 3

rescale3 <- c('Subtype_specificity_NSCLC', 'Subtype_specificity_SCLC', 'Subtype_specificity_GBM', 'Subtype_specificity_MBL')
data[, rescale3] <- (data[, rescale3] /3.069652 ) # Largest value in the dataset is 3.069652, so divide by 3.069652 to get values between 0 and 1

# Identify the columns you want to set values to 0 in
columns_to_set_to_0 <- c(
  "Enrichment_drivers_NSCLC", "Enrichment_drivers_SCLC", "Enrichment_drivers_GBM", "Enrichment_drivers_MBL",
  "Enrichment_mutations_NSCLC", "Enrichment_mutations_SCLC", "Enrichment_mutations_GBM", "Enrichment_mutations_MBL"
)

# Use ifelse to set values less than 0 to 0 in the specified columns
for (col in columns_to_set_to_0) {
  data[, col] <- ifelse(data[, col] < 0, 0, data[, col])
}


palettes2 <- tribble(
  ~palette,             ~colours,          ~limits,
  "overall",            grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Greys")[-1]))(101), c(0, 1),
  "palette1",           grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues") %>% c("#011636"))(1001), c(0.7, 1),
  "palette2",           grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Reds")[-8:-9])(1001), c(0.85, 1),
  "palette3",           grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Greens")[-8:-9])(1001), c(0, 0.25)
)


column_groups <- tribble( # tribble_start
  ~Category,  ~group,         ~palette,
  "Network Construction",  "Network Construction",      "overall",
  "Hub\nspecificity",  "Hub specificity",      "benchmark",
  "IntOGen\nenrichment",  "Enrichment known drivers",       "benchmark2",
  
  "CMP\nenrichment",  "Enrichment mutations",       "benchmark2",
  "Omics\ncorrelation",  "Omics correlation",       "scaling",
) # tribble_end


legends <- list(list(title = "Specificity",
                     palette = "benchmark",
                     geom = "circle",
                     labels=c("0.0","", "0.5","", "1")),
                list(title = "1 - p-value",
                     palette = "benchmark2",
                     geom = "circle",
                     labels = c("0.5","", "0.75","", "1")),
                list(title = "PCC with omics",
                     palette = "scaling",
                     geom = "circle",
                     labels = c("0.05","", "0.5","", "1"))
                )

column_info <- tribble(
  ~id,     ~group,         ~name,                      ~geom,        ~palette,    ~options,
  "Tools",    "",             "",                         "text",       NA,          list(hjust = 0, width = 3.5),
  "Aggregate_network",    "Network Construction",      "Aggregate network",         "text",       NA,          list(hjust = 0.5, width = 1.5),
  "Background_network",    "Network Construction",      "Background network",         "text",       NA,          list(hjust = 0.5, width = 1.5),
  "Scalable",    "Network Construction",      "Scalable",         "text",       NA,          list(hjust = 0.5, width = 1.5),
  "Implementation",    "Network Construction",      "Implementation",         "text",       NA,          list(hjust = 0.5, width = 3),
  "Allows_advanced_GRN_inference",    "Network Construction",      "Other GRN tools",         "text",       NA,          list(hjust = 0.5, width = 2),
  "Edge_bias",    "Network Construction",      "Edge weight bias",         "text",       NA,          list(hjust = 0.5, width = 2),
  
  
  "Subtype_specificity_NSCLC",    "Hub specificity",      "NSCLC",         "circle",       "benchmark",          lst(width = 1.2),
  "Subtype_specificity_SCLC",    "Hub specificity",      "SCLC",         "circle",       "benchmark",          lst(width = 1.2),
  "Subtype_specificity_GBM",    "Hub specificity",      "GBM",         "circle",       "benchmark",          lst(width = 1.2),
  "Subtype_specificity_MBL",    "Hub specificity",      "MBL",         "circle",       "benchmark",          lst(width = 1.2),
  
  
  "Enrichment_drivers_NSCLC",    "Enrichment known drivers",      "NSCLC",         "circle",       "benchmark2",          lst(width = 1.2),
  "Enrichment_drivers_SCLC",    "Enrichment known drivers",      "SCLC",         "circle",       "benchmark2",          lst(width = 1.2),
  "Enrichment_drivers_GBM",    "Enrichment known drivers",      "GBM",         "circle",       "benchmark2",          lst(width = 1.2),
  "Enrichment_drivers_MBL",    "Enrichment known drivers",      "MBL",         "circle",       "benchmark2",          lst(width = 1.2),
  
  "Enrichment_mutations_NSCLC",    "Enrichment mutations",      "NSCLC",         "circle",       "benchmark2",          lst(width = 1.2),
  "Enrichment_mutations_SCLC",    "Enrichment mutations",      "SCLC",         "circle",       "benchmark2",          lst(width = 1.2),
  "Enrichment_mutations_GBM",    "Enrichment mutations",      "GBM",         "circle",       "benchmark2",          lst(width = 1.2),
  "Enrichment_mutations_MBL",    "Enrichment mutations",      "MBL",         "circle",       "benchmark2",          lst(width = 1.2),
  
  "Correlation_lung_proteomics",    "Omics correlation",      "Lung proteomics",         "circle",       "scaling",          lst(width = 1.2),
  "Correlation_lung_CNV",    "Omics correlation",      "Lung CNV",         "circle",       "scaling",          lst(width = 1.2),
  "Correlation_brain_proteomics",    "Omics correlation",      "Brain proteomics",         "circle",       "scaling",          lst(width = 1.2),
  "Correlation_brain_CNV",    "Omics correlation",      "Brain CNV",         "circle",       "scaling",          lst(width = 1.2),
)

funky_heatmap(data, column_info = column_info, column_groups = column_groups, expand = list(xmax = 4), palettes = palettes2, legends = legends, scale = FALSE)

# funky_heatmap(data, column_info = column_info, column_groups = column_groups, expand = list(xmax = 4), palettes = palettes3, scale = FALSE)
g <- funky_heatmap(data, column_info = column_info, column_groups = column_groups, expand = list(xmax = 4), palettes = palettes2, legends = legends, scale = FALSE)
#funky_heatmap(data, column_info = column_info, column_groups = column_groups, expand = list(xmax = 4), palettes = palettes2, scale = TRUE)
                
#g <- funky_heatmap(data, column_info = column_info, column_groups = column_groups, expand = list(xmax = 4), palettes = palettes2)

ggsave("/home/boris/Documents/PhD/single_sample/results/Rebuttal/SWEET/Funkyheatmap2.pdf", g, device = cairo_pdf, width = g$width, height = g$height)
