#### Running DGE for GSE193123 and generating plots of the data

# Installing necessary packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("ggplot2")
install.packages("pheatmap")
BiocManager::install("EnhancedVolcano")
BiocManager::install("tximport")
BiocManager::install("sva")
install.packages("reshape2")
BiocManager::install("AnnotationDbi")
BiocManager::install("EnsDb.Hsapiens.v75")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("goseq")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("biomaRt")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install(c("GenomicFeatures", "IRanges"))

# Loading necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(tximport)
library(RColorBrewer)
library(reshape2)
library(sva)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v75)
library(org.Hs.eg.db)
library(goseq)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(IRanges)
library(ggplot2)


# Set working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/RNAseq data")

### Setting up data
# Read in count matrix 
count_data <- read.table("gene_counts.txt", header = TRUE, row.names = 1, sep = "\t")

# Check the structure of the dataset
class(count_data)

# Remove the first five columns that contain non-sample data
count_data <- count_data[, 6:25]

### Defining sample information
sample_info <- data.frame(
  sample = colnames(count_data),  
  condition = c("control", "control", "control", "control", "control", "control", "control", 
                "control", "control", "control", "PCOS", "PCOS", "PCOS", "PCOS", "PCOS","PCOS",
                "PCOS", "PCOS", "PCOS", "PCOS")
)

# Convert condition to a factor
sample_info$condition <- factor(sample_info$condition, levels = c("control", "PCOS"))

# Set row names as sample names
rownames(sample_info) <- sample_info$sample

# Remove sample column since row names represent samples
sample_info$sample <- NULL  



### Checking data and saving file
# Final check: Does count_data still have numerical values?
print(head(count_data))
print(sample_info)

# Ensure row names in sample_info match colnames in count_data
all(rownames(sample_info) %in% colnames(count_data))  # TRUE

# Save cleaned data
write.csv(count_data, "cleaned_gene_counts.csv", row.names = TRUE)
write.csv(sample_info, "sample_info.csv", row.names = TRUE)




# Running differential gene expression analysis and saving results
# Create DESeq2 dataset to work with
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ condition)

# Running the DGE
dds <- DESeq(dds)

# Getting results
res <- results(dds, contrast = c("condition", "PCOS", "control"))

# Convert DESeq2 results to a dataframe
res_df <- as.data.frame(res)  

# Preserve gene IDs
res_df$GeneID <- rownames(res_df) 

# Remove Gene_ID column
rownames(res_df) <- NULL  

# Saving the results in .csv
write.csv(res_df, file = "results_unfiltered.csv", row.names = FALSE)




### Filtering results to get top significant genes 
# Apply filtering based on adjusted p-value and log2FoldChange
sig_res <- res_df[which(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.263), ]



### Mapping Entrez IDs and saving file
# Mapping Entrez IDs for significant results (sig_res)
sig_res$entrez <- mapIds(org.Hs.eg.db, 
                         keys = rownames(sig_res), 
                         column = "ENTREZID", 
                         keytype = "ENSEMBL", 
                         multiVals = "first")

# Removing rows where Entrez ID is missing 
sig_res <- sig_res[!is.na(sig_res$entrez), ]

# Sorting by adjusted p-value 
sig_res <- sig_res[order(sig_res$padj), ]

# So now, `sig_res` contains filtered, significant results with Entrez IDs

# Saving the results in .csv
write.csv(sig_res, file = "filtered_significant_genes.csv", row.names = FALSE)





#### Transformation and normalization 
# Applying Variance Stabilizing Transformation (VST) for normalized counts
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# Extracting the normalized counts from VST
norm_counts <- assay(vsd)

# Select the top 50 DEGs based on adjusted p-value 
top50_genes <- rownames(sig_res)[1:50]  
top50_counts <- norm_counts[top50_genes, ]




### Visualization of DEG
# 1. Creating Heatmap
pdf("heatmap_top50_DEGs.png", width = 800, height = 600)
pheatmap(top50_counts, 
         cluster_rows = TRUE, cluster_cols = TRUE, 
         scale = "row",
         show_rownames = TRUE, annotation_col = sample_info,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         main = "Heatmap of Top 50 DEGs")
dev.off()



# 2. Creating volcano plot
sig_res <- fread("filtered_significant_genes.csv")

# Create a data frame with results for volcano plot
volcano_data <- as.data.frame(res)
volcano_data$significant <- ifelse(volcano_data$padj < 0.05 & abs(volcano_data$log2FoldChange) >= 1, "Yes", "No")

# Create the volcano plot
png("volcano_plot.png", width = 800, height = 600)
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("No" = "gray", "Yes" = "red")) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Add threshold for p-value
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")  # Add threshold for fold change
dev.off()





