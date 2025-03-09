###### Pre-processing for Methyseq data E-MTAB-3777

install.packages("tidyverse")


library(readr)
library(data.table)
library(dplyr)
library(tidyverse)

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
#all(rownames(meth_a) == rownames(meth_b))  #Returned true

# Check for missing values
# sum(is.na(meth_a)) 
# sum(is.na(meth_b))

# Check for duplicate sample names
#colnames(meth_b)[!(colnames(meth_b) %in% colnames(meth_a))]  ##"ILMNID", "NAME", and "V114"  looked pretty different              

# Check if "ILMNID", "NAME", and "V114" contain redundant information
# head(meth_b[, c("ILMNID", "NAME", "V114")]) 

# Check the differences in column names between the two samples
# setdiff(colnames(meth_a), colnames(meth_b))   # column v112 is unique
# setdiff(colnames(meth_b), colnames(meth_a))   # ILMNID and Names are redundant but column V114 is unique

# Check content of V112 and V114
# meth_a$V112  # contains N/A values
# meth_b$V114  # contains N/A

# Check if all the values in v112 and v114 are N/A
# all(is.na(meth_a$V112))  # TRUE
# all(is.na(meth_b$V114))  # TRUE

# Remove column V112 from meth_a
meth_a <- meth_a %>% select(-V112)

# Remove the 3 redundant columns from meth_b
meth_b <- meth_b %>% select(-c(ILMNID, NAME, V114))


# 3. Merge the dataframes by TargetID, probeID_A, and probeID_B
merged_data <- inner_join(meth_a, meth_b, by = c("TargetID", "ProbeID_A", "ProbeID_B"))

# Set the TargetID (CpG sites) as row names 
rownames(merged_data) <- meth_a$TargetID

# Remove the TargetID 
merged_data$TargetID <- NULL

# Save merged_data as a file
write.csv(merged_data, "merged_meth.csv", row.names = TRUE)


# 4. Getting the CpG x sample matrix
# Extract sample IDs
column_names <- colnames(merged_data)
sample_ids <- unique(str_extract(column_names, "^\\d+"))
sample_ids <- sample_ids[!is.na(sample_ids)]

# Extract AVG_Beta Columns and Tidy
tidy_beta <- merged_data %>%
  mutate(cpg_site = rownames(merged_data)) %>%
  pivot_longer(cols = -cpg_site, names_to = "variable", values_to = "AVG_Beta") %>%
  filter(str_detect(variable, "AVG_Beta")) %>%
  mutate(sample_id = str_extract(variable, "^\\d+")) %>%
  select(cpg_site, sample_id, AVG_Beta)

# Create a wide dataframe.
wide_beta <- tidy_beta %>%
  pivot_wider(names_from = sample_id, values_from = AVG_Beta)

# Convert to dataframe
wide_beta <- as.data.frame(wide_beta)  

# Reset the TargetID (CpG sites) as row names 
rownames(wide_beta) <- wide_beta$cpg_site

# Remove the TargetID 
wide_beta$cpg_site <- NULL

# Save sample-by-CpG matrix for analyses
write.csv(wide_beta, "meth_matrix.csv", row.names = TRUE)
