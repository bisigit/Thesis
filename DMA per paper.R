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
library(pheatmap)
library(tidyr)
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

# save Results
write.csv(beta_filtered, "beta_vals_for_DMA.csv", row.names = TRUE) 


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
metadata$`Clinical History` <- factor(metadata$`Clinical History`)

# Define design matrix
design <- model.matrix(~ group + Age + BMI + `Clinical History`, data = metadata)

# Set rownames to sample IDs
rownames(design) <- colnames(m_imputed)


# Fit the model
fit <- lmFit(m_imputed, design)
fit <- eBayes(fit)

# Full table, adjusted using Benjamini-Hochberg
dmp_results <- topTable(fit, coef = "groupPCOS", number = Inf, adjust = "BH")    ## Significant results found


# Filter results based on padj and log2foldchange parameters
res_sig <- 
  dmp_results[which(dmp_results$adj.P.Val < 0.05 & abs(dmp_results$logFC) > 0.263),]

# Add a column that indicates up and downregulated genes
res_sig$direction <- ifelse(res_sig$logFC > 0, "Up", "Down")
table(res_sig$direction)   # Down 18,  Up 2

# save Results
write.csv(res_sig, "dmp_results latest.csv", row.names = TRUE)


# Getting list of DMR
dmgs <- row.names(res_sig)


# save Results
write.csv(dmgs, "dmp_sig_list.csv", row.names = TRUE)




################## Clustered Heatmap of DEGs #######################
# Extract significant probes
sig_probes <- rownames(res_sig)

# Subset the methylation matrix
heatmap_data <- beta_filtered[sig_probes, ]  # beta values (0–1 scale)


# Check probe names
intersect(sig_probes, rownames(beta_filtered))


# Add annotations
annotation_col <- data.frame(Group = metadata$Disease)
rownames(annotation_col) <- colnames(heatmap_data)


# Plot and save heatmap
png("methylation_heatmap.png", width = 1000, height = 1200, res = 150)
pheatmap(
  heatmap_data,                          # gene-level beta or M-values
  color = colorRampPalette(c("blue", "white", "red"))(50),
  annotation_col = annotation_col,     # e.g., PCOS/control
  scale = "row",                       # normalize across samples
  cluster_cols = FALSE,
  show_rownames = FALSE,
  main = "DMGs: Differential Methylation"
)
dev.off()






###################### Annotation for GSEA ##############################
# Annotate cpgs
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
res_sig$TargetID <- rownames(res_sig)
annotated_dmps <- merge(res_sig, ann, by.x = "TargetID", by.y = "Name")


# 1. Identify Genes with Significant Methylation Changes
gene_list <- unique(unlist(strsplit(annotated_dmps$UCSC_RefGene_Name, split = ";|,")))


# 2. Subset by genomic region
promoter_dmps <- annotated_dmps[grep("TSS1500|TSS200", annotated_dmps$UCSC_RefGene_Group), ]

# 3. Check direction of change
annotated_dmps$methylation_direction <- ifelse(annotated_dmps$logFC > 0, "Hypermethylated", "Hypomethylated")
table(annotated_dmps$methylation_direction)

# 4. Prep for enrichment
sig_cpgs <- annotated_dmps$TargetID


# Save annotated results
write.csv(annotated_dmps, "annotated_dmp_results.csv", row.names = TRUE)

















########################## Exploratory Analysis ###################################
# 1. Plot Heatmap of top 500 most variable genes
# Reconvert to beta values
beta_vis <- 2^m_imputed / (2^m_imputed + 1)

# Compute row variance
row_var <- apply(beta_vis, 1, var)

# Select top 500 most variable CpGs
top500_cpgs <- names(sort(row_var, decreasing = TRUE))[1:500]

# Subset the B-values matrix
top500_matrix <- beta_vis[top500_cpgs, ]


# 2. Plot Heatmap to see which genes are variable
# Define annotation for columns (samples)
annotation_col <- data.frame(
  Group = metadata$Disease
)
rownames(annotation_col) <- metadata$`Source No`

# Heatmap
heatmap_explore <- pheatmap(
  top500_matrix,
  scale = "row",  # standardize CpGs
  show_rownames = FALSE,
  annotation_col = annotation_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete"
)


# 3. Boxplots of Top DMPs (Using β-values)
# Extract top 6 DMPs by unadjusted p-value (or adjusted if any are significant)
top_cpgs <- rownames(dmp_results)[1:10]

# Loop to plot each CpG
for (cpg in top_cpgs) {
  plot_data <- data.frame(
    Beta = beta_vis[cpg, ],
    Group = metadata$Disease,
    Sample = colnames(beta_vis)
  )
  
  p <- ggplot(plot_data, aes(x = Group, y = Beta, fill = Group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, shape = 21, size = 2) +
    labs(title = paste("β-value for", cpg), y = "Beta (Methylation %)") +
    theme_minimal()
  
  print(p)  # Or save using ggsave
}
ggsave(paste0("beta_plot_", cpg, ".png"), p, width = 6, height = 4, dpi = 300)





# Print top 10 DMRs to have them
print(top_cpgs)
# [1] "cg05999287" "cg08567382" "cg25504273" "cg24120153" "cg26890920" "cg26807939" "cg17977480" "cg12533748"
# [9] "cg21736038" "cg21184320"



























