#### Down-sampling GSE84958 ####

library(ggplot2)
library(readr)
library(data.table)
library(dplyr)

# Set working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/RNAseq data")

# 1. Download Metadata from SRA Run Selector
# Go to GEO or SRA Run Selector and download the metadata as a CSV file (SraRunTable.csv)
# Separate Run file into control samples and PCOS samples

# Load control and PCOS metadata separately
controls <- read.csv("controls.csv", header=TRUE)
pcos <- read.csv("pcos.csv", header=TRUE)

# Apply filtering conditions
filtered_controls <- controls %>%
  filter(AvgSpotLen > 75, Platform == "ILLUMINA", LibraryLayout == "PAIRED") %>%
  sample_n(10)  # Randomly select 10 controls

filtered_pcos <- pcos %>%
  filter(AvgSpotLen > 75, Platform == "ILLUMINA", LibraryLayout == "PAIRED") %>%
  sample_n(10)  # Randomly select 10 PCOS

# Combine selected samples
selected_samples <- bind_rows(filtered_controls, filtered_pcos)

# Save selected SRA Run IDs for downloading
write.table(selected_samples$Run, "selected_samples.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# View final selection
selected_samples
