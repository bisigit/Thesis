########### DMR Analysis ###############

# Load Libraries
library(readr)
library(data.table)
library(dplyr)
library(tidyverse)
library(minfi)
library(limma)
library(impute)
library(mice)
library(sva)
library(ggplot2)
library(preprocessCore)
library(oligo)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Set working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/METHarray data")

# Load matrix
wide_beta <- fread("meth_matrix.csv", header = TRUE)
rownames(wide_beta) <- wide_beta$V1
wide_beta$V1 <- NULL

# Load metadata
metadata <- fread("E-MTAB-3777.sdrf.txt", header = TRUE)

# Rename disease states
metadata$`Characteristics[disease]` <- recode(metadata$`Characteristics[disease]`, 
                                    "normal" = "control", 
                                    "polycystic ovary syndrome" = "PCOS")

# Match sample order in methylation matrix with metadata
metadata <- metadata[match(colnames(wide_beta), metadata$`Source Name`), , drop = FALSE]




# A. QC
# 1. Remove probes on the y chromosome
# Load annotation data
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Extract probe names that are NOT on the Y chromosome
keep_probes <- annotation$Name[annotation$chr != "chrY"]
beta_filtered <- wide_beta[rownames(wide_beta) %in% keep_probes, ]

# Ensure row names are retained
rownames(beta_filtered) <- rownames(wide_beta)[rownames(wide_beta) %in% keep_probes]

# B. Convert β-Values to M-Values
m_values <- log2(beta_filtered / (1 - beta_filtered))

# Apply quantile normalization
m_norm <- normalize.quantiles(as.matrix(m_values))
rownames(m_norm) <- rownames(beta_filtered)
colnames(m_norm) <- colnames(beta_filtered)
m_norm <- as.data.frame(m_norm)


# Handle NA values
# 1. Remove rows with more that 50% missingness
# Set threshold for missing values (50%)
threshold <- 0.5 * ncol(m_norm)

# Identify rows to keep
rows_to_keep <- rowSums(is.na(m_norm)) <= threshold

# Subset while preserving row names
m_norm_filtered <- m_norm[rows_to_keep, , drop = FALSE]

# 2. Apply KNN imputation
m_norm_imputed <- impute.knn(as.matrix(m_norm_filtered))$data




# C. SVA to correct for age and BMI differences
# Define group, age, and BMI
group <- factor(metadata$`Characteristics[disease]`, levels = c("control", "PCOS"))
age <- metadata$`Characteristics[age]`
bmi <- metadata$`Characteristics[bmi]`

# Define model matrices
mod <- model.matrix(~ group + age + bmi)
mod0 <- model.matrix(~ age + bmi)

# Estimate surrogate variables
n.sv <- num.sv(m_norm_imputed, mod, method = "be")

# Run SVA
sva_results <- sva(m_norm_imputed, mod, mod0, n.sv = n.sv)

# Extract surrogate variables
sv_matrix <- sva_results$sv

# Update design matrix to include surrogate variables
design <- cbind(mod, sv_matrix)



# D. DMA
# Fit linear model
fit <- lmFit(m_norm_imputed, design)
fit <- eBayes(fit)

# Extract results for PCOS vs. Control
dmp_results <- topTable(fit, coef = "groupPCOS", number = Inf, adjust = "BH")

# Filter significant sites (adj. p-value < 0.05)
dmp_significant <- dmp_results[dmp_results$adj.P.Val < 0.05, ]   # No significant regions

write.csv(dmp_results, "dmp.csv", row.names = TRUE)
