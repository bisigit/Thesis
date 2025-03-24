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
library(DMRcate)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Set working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/METHarray data")

# 1. Load matrix
wide_beta <- fread("filtered_meth_matrix.csv", header = TRUE)
rownames(wide_beta) <- wide_beta$V1
wide_beta$V1 <- NULL

# Rename columns with sample IDs
colnames(wide_beta) <- sub("\\..*", "", colnames(wide_beta))
beta_matrix <- as.data.frame(wide_beta)
row.names(beta_matrix) <- row.names(wide_beta)

# 2. Subsetting data used in paper
# List of samples in paper
samples_used <- c(
  "90857", "90966", "91121", "91140", "91364", 
  "91410", "91414", "91452", "91512", "91520",  # PCOS
  "90802", "90590", "90860", "91032", "91066", 
  "91225", "91228", "91301", "91443", "91517"   # Controls
)

# Subset data
samples_in_data <- intersect(samples_used, colnames(beta_matrix))  
beta_subset <- beta_matrix[, samples_in_data]

#  Load metadata
metadata <- fread("metadata.csv", header = TRUE)

# Subset and re-order metadata to match beta matrix
metadata <- metadata[match(colnames(beta_subset), metadata$`Source No`), ]

# Check that metadata and beta matrix are perfectly aligned
stopifnot(all(metadata$`Source No` == colnames(beta_subset)))



# 3. Remove probes on the y chromosome
# Load annotation data
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Extract probe names that are NOT on the Y chromosome
keep_probes <- annotation$Name[annotation$chr != "chrY"]
beta_filtered <- beta_subset[rownames(beta_subset) %in% keep_probes, ]

# Ensure row names are retained
rownames(beta_filtered) <- rownames(beta_filtered)[rownames(beta_filtered) %in% keep_probes]


# 4. Convert β-Values to M-Values
offset <- 1e-6
m_values <- log2((beta_filtered + offset) / (1 - beta_filtered + offset))


# 5. Handle Missing values
# Set threshold for missing values (50%)
threshold <- 0.5 * ncol(m_values)

# Identify rows to keep
rows_to_keep <- rowSums(is.na(m_values)) <= threshold

# Subset while preserving row names
m_values_filtered <- m_values[rows_to_keep, , drop = FALSE]



# 6. Apply KNN imputation
m_imputed <- impute.knn(as.matrix(m_values_filtered))$data



# 7. DMA
# Define group variable
group <- factor(metadata$Disease, levels = c("control", "PCOS"))

# Define deisign matrix
design <- model.matrix(~ group) 

# Set rownames to sample IDs
rownames(design) <- colnames(m_imputed)


# Fit the model
fit <- lmFit(m_imputed, design)
fit <- eBayes(fit)

# Full table, adjusted using Benjamini-Hochberg
dmp_results <- topTable(fit, coef = "groupPCOS", number = Inf, adjust = "BH")    ## No significant results




# 8. DMR Analysis
# Ensure design rownames match colnames of matrix

# Ensure alignment
stopifnot(identical(colnames(m_imputed), rownames(design)))

# Run DMR
annotated <- cpg.annotate(
  object = m_imputed, 
  datatype = "array", 
  what = "M", 
  analysis.type = "differential", 
  design = design, 
  coef = "groupPCOS",
  arraytype = "450K"
)

# Your contrast returned no individually significant probes. 
# Try increasing the fdr. Alternatively, set pcutoff manually in dmrcate() 
# to return DMRs, but be warned there is an increased risk of Type I errors.

# Increase the pcutoff
dmrcoutput <- dmrcate(annotated, lambda = 1000, C = 2, pcutoff = 0.2)
results.ranges <- extractRanges(dmrcoutput)

# Get results table
dmr_table <- as.data.frame(results.ranges)
View(dmr_table)

# Filter by mean difference
top_dmrs <- dmr_table[order(-abs(dmr_table$meandiff)), ]

# Extract overlapping genes
promoters <- as.character(results.ranges$overlapping.promoters)
genes <- unique(unlist(strsplit(promoters, ";")))
head(genes)   #didn't work

