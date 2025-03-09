#### DMA for Methyseq data E-MTAB-3777



library(data.table)
library(minfi)
library(limma)
library(ggplot2)
library(mice)
library(EnhancedVolcano)
library(sva)



# Setting work directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/METHarray data")

# 1. Load data
completed_data <- fread("matrix_for_dma.csv", header = TRUE)

# Set the Cpg_id column as row names
row.names(completed_data) <- completed_data$V1
completed_data$V1 <- NULL

# Load metadata file
metadata <- fread("E-MTAB-3777.sdrf.txt", header = TRUE)

# Extract sample names
sample_names <- colnames(completed_data)

# 2. PCA
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


# 3. Remove outlier
# Remove the outlier column
completed_data$`polycystic ovary syndrome.2` <- NULL

# Remove the row in sample_names corresponding to the outlier
sample_names <- sample_names[sample_names != "polycystic ovary syndrome.2"]


# 4. SVA analysis
# Making sample_names a dataframe
# Create a metadata dataframe
sample_names <- data.frame(
  SampleName = sample_names,
  DiseaseState = ifelse(grepl("polycystic", sample_names), "PCOS", "Control") 
)


# Define the model matrix (biological variable of interest)
mod <- model.matrix(~ DiseaseState, data = sample_names)

# Define the null model 
mod0 <- model.matrix(~ 1, data = (sample_names))

# Estimate the number of surrogate variables
num.sv <- num.sv(as.matrix(completed_data), mod, method = "be")

# Run SVA to obtain surrogate variables
sva_results <- sva(as.matrix(completed_data), mod, mod0, n.sv = num.sv)

# Add surrogate variables to the design matrix
mod_sv <- cbind(mod, sva_results$sv)

# Fit the linear model with surrogate variables included
fit <- lmFit(completed_data, mod_sv)
fit <- eBayes(fit)

# Extract top differentially methylated positions (DMPs)
dmps <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "P")


# Save SVA-corrected data
write.csv(completed_data, "sva_corrected_matrix.csv", row.names = TRUE)


# 5. Re-run PCA after outlier removal and SVA correction
# PCA analysis after SVA correction
pca_result_2 <- prcomp(t(completed_data), center = TRUE, scale. = TRUE)

# Convert PCA results to a dataframe
pca_data <- data.frame(pca_result_2$x)
pca_data$SampleID <- colnames(completed_data)

# Scatter plot of the first two principal components
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, label = SampleID)) +
  geom_point(size = 5) +
  geom_text(vjust = -1, size = 3) +  # Label points with sample names
  labs(title = "PCA of Methylation Data (Post-SVA)", x = "PC1", y = "PC2") +
  theme_minimal()

# Save PCA plot
ggsave("pca_plot_2.png", plot = pca_plot, width = 18, height = 8)
ggsave("pca_plot_2.pdf", plot = pca_plot, width = 18, height = 8)  # Still no cluster but outlier gone


# 6. DMA with Limma
# Define sample groups 
group <- factor(sample_names$DiseaseState)
design <- model.matrix(~ group)

# Fit linear model
fit <- lmFit(completed_data, design)
fit <- eBayes(fit)

# Get differentially methylated positions (DMPs)
dmps <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "P")

# Apply thresholds for FDR and log2 fold-change
dmps_filtered <- dmps[dmps$adj.P.Val < 0.05 & abs(dmps$logFC) > 0.263, ]  #NO DMRs 


