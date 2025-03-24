############### Using Signals #######################

library(readr)
library(data.table)
library(minfi)  
library(dplyr)  
library(wateRmelon)
library(preprocessCore)

# Set working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/METHarray data")

# 1. Manually download dataset file the files section in ArrayExpress

# 2. Read the processed methylation data
meth_a <- fread("PCOS_SubQ_Adipose_Methylation_450k_a.txt", header = TRUE)
meth_b <- fread("PCOS_SubQ_Adipose_Methylation_450k_b.txt", header = TRUE)
metadata <- fread("E-MTAB-3777.sdrf.txt", header = TRUE)

# Convert to data frame
meth_a <- as.data.frame(meth_a)
meth_b <- as.data.frame(meth_b)

# Check if CpG IDs match between files
all(rownames(meth_a) == rownames(meth_b))  

# Check for duplicate sample names
colnames(meth_b)[!(colnames(meth_b) %in% colnames(meth_a))]  ##"ILMNID", "NAME", and "V114"  looked pretty different from the other columns

# Check if "ILMNID", "NAME", and "V114" contain redundant information
head(meth_b[, c("ILMNID", "NAME", "V114")]) #Yes

# Check the differences in column names between the two samples
setdiff(colnames(meth_a), colnames(meth_b))   # column v112 is unique
setdiff(colnames(meth_b), colnames(meth_a))   # ILMNID and Names are redundant but column V114 is unique

# Check content of V112 and V114
meth_a$V112  # contains N/A values
meth_b$V114  # contains N/A

# Check if all the values in v112 and v114 are N/A
all(is.na(meth_a$V112))  # TRUE
all(is.na(meth_b$V114))  # TRUE

# Remove column V112 from meth_a
meth_a <- meth_a %>% select(-V112)

# Remove the 3 redundant columns from meth_b
meth_b <- meth_b %>% select(-c(ILMNID, NAME, V114))

# 3. Merge the dataframes by TargetID, probeID_A, and probeID_B
merged_data <- inner_join(meth_a, meth_b, by = c("TargetID", "ProbeID_A", "ProbeID_B"))




######################### Filtering #######################
# 1. Extract all relevant probes 
# Extract detection p-values
detection_pvals <- merged_data[, grep("Detection Pval", colnames(merged_data))]

# Extract bead counts
bead_counts_A <- merged_data[, grep("Avg_NBEADS_A", colnames(merged_data))]
bead_counts_B <- merged_data[, grep("Avg_NBEADS_B", colnames(merged_data))]

# Extract bead standard errors
bead_stderr_A <- merged_data[, grep("BEAD_STDERR_A", colnames(merged_data))]
bead_stderr_B <- merged_data[, grep("BEAD_STDERR_B", colnames(merged_data))]


# 2. Filter Out Low-Quality Probes
# Remove probes with high detection values
keep_probes <- rowMeans(detection_pvals < 0.01) > 0.75
meth_filtered <- merged_data[keep_probes, ]

# Remove probes with low bead count
keep_beads <- rowMeans(bead_counts_A >= 3) > 0.75 & rowMeans(bead_counts_B >= 3) > 0.75

# Apply filtering while ensuring unique TargetID
meth_filtered <- meth_filtered[keep_beads, ]

# Remove potential duplicate TargetIDs
meth_filtered <- meth_filtered[!duplicated(meth_filtered$TargetID), ]

# 3. Remove Probes with High Bead Standard Error
# Find the 90th percentile 
stderr_threshold_A <- quantile(bead_stderr_A, 0.99, na.rm = TRUE)  
stderr_threshold_B <- quantile(bead_stderr_B, 0.99, na.rm = TRUE)

# Count probes above the threshold
# Check if any value in each row exceeds the threshold
exceed_threshold <- rowSums(bead_stderr_A > stderr_threshold_A | bead_stderr_B > stderr_threshold_B, na.rm = TRUE) > 0

# Count how many probes exceed the threshold
sum(exceed_threshold)  # [1] 114546

# Filter out the probes exceeding the threshold
meth_filtered <- meth_filtered[!exceed_threshold, ]

# Remove potential duplicate TargetIDs
meth_filtered <- meth_filtered[!duplicated(meth_filtered$TargetID), ]
meth_filtered <- meth_filtered[!is.na(meth_filtered$TargetID), ]



########################### B values from raw intensities ###################
# 1. Obtain signal intensities
signal_unmeth <- meth_filtered[, grep("Signal_A", colnames(meth_filtered))]
signal_meth   <- meth_filtered[, grep("Signal_B", colnames(meth_filtered))]

# 2. Normalize raw intensities
# Quantile normalize each channel separately
signal_unmeth_norm <- normalize.quantiles(as.matrix(signal_unmeth))
signal_meth_norm   <- normalize.quantiles(as.matrix(signal_meth))

# Restore column names
# Restore row/column names
colnames(signal_unmeth_norm) <- colnames(signal_unmeth)
colnames(signal_meth_norm)   <- colnames(signal_meth)


# 3. Calculate beta values
beta_matrix <- signal_meth_norm / (signal_meth_norm + signal_unmeth_norm + 100)
rownames(beta_matrix) <- meth_filtered$TargetID

# 5. Save sample-by-CpG matrix for analyses
write.csv(beta_matrix, "Filtered_meth_matrix.csv", row.names = TRUE)

# 4. Convert to M-values
#offset <- 1e-6
#m_values <- log2((beta_matrix + offset) / (1 - beta_matrix + offset))

# 5. Save sample-by-CpG matrix for analyses
#write.csv(beta_matrix, "Filtered_meth_matrix.csv", row.names = TRUE)

