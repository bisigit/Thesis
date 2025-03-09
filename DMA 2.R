#### DMA for Methyseq data E-MTAB-3777

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi")
install.packages("mice")
BiocManager::install("EnhancedVolcano")

library(data.table)
library(minfi)
library(limma)
library(ggplot2)
library(mice)
library(EnhancedVolcano)
library(sva)



# Setting work directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/METHarray data")

# 1. Import a CSV file into R
meth_matrix <- fread("meth_matrix.csv", header = TRUE)

# Load metadata file
metadata <- fread("E-MTAB-3777.sdrf.txt", header = TRUE)

# Set the Cpg_id column as row names
row.names(meth_matrix) <- meth_matrix$V1
meth_matrix$V1 <- NULL

# Transform the beta values to M-values (log transformation)
m_value_matrix <- log2(meth_matrix / (1 - meth_matrix))

# Assign CpG IDs from the original beta_values matrix to the m_values matrix
rownames(m_value_matrix) <- rownames(meth_matrix)

# Checking before normalization 
boxplot(m_value_matrix, main = "M-values Before Normalization", outline = FALSE)

# 2. Apply quantile normalization (paper states that limma package was used for data analysis)
normalized_m_values <- normalizeQuantiles(m_value_matrix)

# Checking data after normalization
boxplot(normalized_m_values, main = "M-values After Normalization", outline = FALSE)

# Check for NA values
sum(is.na(meth_matrix))

# Calculate the percentage of NA values
na_percentage <- sum(is.na(normalized_m_values)) / (nrow(normalized_m_values) * ncol(normalized_m_values)) * 100
print(paste("Percentage of NA values:", na_percentage, "%"))    # [1] "Percentage of NA values: 0.14476763383219 %"

# 3. Carryout imputation to fix missing values
# Imputation using mice
imputed_data <- mice(normalized_m_values, m = 1, method = "pmm", seed = 123)
completed_data <- complete(imputed_data, 1)

# 4. Rename column names to reflect the grouping
#sample_names <- data.frame(
#  sample_ID = c("90242", "90590", "90666", "90802", "90857", "90860", "90865", "90966", "90971", "91032", "91066", "91121", "91140",
#                "91225", "91228", "91301", "91364", "91410", "91414", "91443", "91452", "91512", "91517", "91520"),
#  condition = c("PCOS", "normal", "PCOS", "normal", "PCOS", "normal", "normal", "PCOS", "PCOS", "normal", "normal", "PCOS", "PCOS", 
#                "normal", "normal", "normal", "PCOS", "PCOS", "PCOS", "normal", "PCOS", "PCOS", "normal", "PCOS"))
sample_names <- metadata[, c("Source Name", "Characteristics[disease]" )]

# Convert SampleID to character 
sample_names$`Source Name` <- as.character(sample_names$`Source Name`)  
colnames(completed_data) <- as.character(colnames(completed_data))  

# Change names of columns
sample_names <- sample_names[match(colnames(completed_data), sample_names$`Source Name`), ]

# Rename the columns of normalized_m_values based on 'Characteristics[disease]'
colnames(completed_data) <- make.unique(sample_names$`Characteristics[disease]`)


# Assign CpG IDs from the original beta_values matrix to the m_values matrix
rownames(completed_data) <- rownames(meth_matrix)

# Save sample-by-CpG matrix for analyses
write.csv(completed_data, "matrix_for_dma.csv", row.names = TRUE)







