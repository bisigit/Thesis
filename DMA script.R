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





# Load data
#completed_data <- fread("matrix_for_dma.csv", header = TRUE)

# Remove first column
#completed_data$V1 <- NULL


# 5. PCA
# PCA is usually performed on scaled data
pca_result <- prcomp(t(completed_data), center = TRUE, scale. = TRUE)

# Summary of PCA to view variance explained by each principal component
summary(pca_result)

# PCA plot (scatter plot of the first two principal components)
pca_data <- data.frame(pca_result$x)
pca_data$SampleID <- colnames(completed_data)

# Scatter plot of the first two principal components
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, label = SampleID)) +
  geom_point(size = 3) +
  geom_text(aes(label = SampleID), vjust = -1, size = 3) +
  labs(title = "PCA of Methylation Data", x = "PC1", y = "PC2") +
  theme_minimal()

ggsave("pca_plot.png", plot = pca_plot, width = 8, height = 6) 
ggsave("pca_plot.pdf", plot = pca_plot, width = 8, height = 6) # No cluster in PCs, 1 outlier. Running SVA analysis


# 6. Remove outlier
# Subset outlier_column 
outlier_column <- "polycystic ovary syndrome.2"

# Remove the outlier column
completed_data <- completed_data[, !colnames(completed_data) %in% outlier_column]

# Remove the row in sample_names corresponding to the outlier
sample_names <- sample_names[-3, ]


# 7. SVA analysis
# Define sample groups 
group <- factor(sample_names$`Characteristics[disease]`)

# Define the model matrix (biological variable of interest)
mod <- model.matrix(~ group)  

# Define the null model 
mod0 <- model.matrix(~ 1, data = sample_names)

# Estimate the number of surrogate variables
num.sv <- num.sv(as.matrix(completed_data), mod, method = "be")

# Run SVA to obtain surrogate variables
sva_results <- sva(as.matrix(completed_data), mod, mod0, n.sv = num.sv)

# Add surrogate variables to the design matrix
sv_matrix <- sva_results$sv %*% t(sva_results$sv)  # Compute batch effect matrix
completed_data <- completed_data - sv_matrix  # Subtract from data

# Save SVA-corrected data
write.csv(completed_data, "sva_corrected_matrix.csv", row.names = TRUE)


# 7. Re-run PCA after outlier removal and SVA correction
# PCA analysis after SVA correction
pca_result_2 <- prcomp(t(completed_data), center = TRUE, scale. = TRUE)

# Convert PCA results to a dataframe
pca_data <- data.frame(pca_result_2$x)
pca_data$SampleID <- colnames(completed_data)

# Scatter plot of the first two principal components
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, label = SampleID)) +
  geom_point(size = 5) +
  geom_text(vjust = -1, size = 3) +  # Label points with sample names
  labs(title = "PCA of Methylation Data", x = "PC1", y = "PC2") +
  theme_minimal()

# Save PCA plot
ggsave("pca_plot_2.png", plot = pca_plot, width = 18, height = 8)
ggsave("pca_plot_2.pdf", plot = pca_plot, width = 18, height = 8)  # Still no cluster and outlier still present



# 8. DMA with Limma
# Define sample groups 
group <- factor(sample_names$`Characteristics[disease]`)
design <- model.matrix(~ group)

# Fit linear model
fit <- lmFit(completed_data, design)
fit <- eBayes(fit)

# Get differentially methylated positions (DMPs)
dmps <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "P")

# Apply thresholds for FDR and log2 fold-change
dmps_filtered <- dmps[dmps$adj.P.Val < 0.05 & abs(dmps$logFC) > 0.263, ]

# View filtered DMPs
head(dmps_filtered) # Empty. No differentially methylated regions












