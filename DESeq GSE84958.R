## DESeq on GSE84958 ##

# Loading relevant libraries
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
library(gplots)



# Set working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/RNAseq data")

# 1. Load and organize data
count_data <- read.table("gene_counts.txt", header = TRUE, row.names = 1, sep = "\t")

# Remove the first five columns that contain non-sample data
count_data <- count_data[, 6:25]

# Filter counts
# Remove genes with less then 10 reads in at least 1 sample
keep <- which(rowSums(count_data > 10) >= 1)
length(keep)
counts <- count_data[keep,]

# 2. Defining sample information
sample_info <- data.frame(
  sample = colnames(counts),  
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

# Ensure row names in sample_info match colnames in count_data
all(rownames(sample_info) %in% colnames(counts))  # TRUE

# Save cleaned data
write.csv(count_data, "cleaned_counts.csv", row.names = TRUE)
write.csv(sample_info, "sample_info.csv", row.names = TRUE)


# 3. Run DESeq
# Create the DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(counts, colData = sample_info, 
                              design = ~condition)


# Run DESeq 
dds <- DESeq(dds)

# Variance stabilization
vsd <- vst(dds, blind = TRUE)


# 4. Running PCA analysis
plotPCA(vsd, intgroup = "condition")


# 5. SVAseq analysis to correct for batch effects
# Defining model matrices
mod <- model.matrix(~ condition, sample_info)
mod0 <- model.matrix(~1, sample_info)

# Run SVAseq
sv <- svaseq(counts(dds), mod, mod0)

# Add surrogate variables to sample_info
sample_info$SV1 <- sv$sv[,1]
sample_info$SV2 <- sv$sv[,2]
sample_info$SV3 <- sv$sv[,3]
sample_info$SV4 <- sv$sv[,4]
sample_info$SV5 <- sv$sv[,5]
sample_info$SV6 <- sv$sv[,6]
sample_info$SV7 <- sv$sv[,7]
sample_info$SV8 <- sv$sv[,8]


# 6. Recreate the dds object with surrogate variables
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + condition)

# 7. Re-Run DESeq
dds <- DESeq(dds)

# Data transformation
vsd <- varianceStabilizingTransformation(dds)


# 8. Visualization
# Get top 500 most variable genes
mads <- rowMads(assay(vsd))
mads <- sort(mads, decreasing = T)
top500 <- mads[1:500]

# Get the matrix of top 500 genes with highest MADs
top500 <- assay(vsd)[which(rownames(assay(vsd)) %in% names(top500)),]
head(top500)

annot <- colData(dds)$condition
levels(annot)
category_colors <- c("control" = "blue", "PCOS" = "red")
mapped_colors <- category_colors[annot]

# Plot a heatmap of top 500 
png("top500_heatmap_clustering.png")
par(mar=c(2,2,2,2))
heatmap.2(top500,
          col=bluered(75),       
          trace="none",         
          dendrogram = "both",
          scale="row",          
          key=TRUE,              
          keysize=1.2,           
          margins=c(6,6),
          labRow = "",
          ColSideColors = mapped_colors
)

# Add a legend for the column annotations
legend("topright",
       legend=c("Normal", "PCOS"), 
       fill=category_colors, 
       title="Diagnosis",
       border=FALSE)
dev.off()



# 7. Detect differentially expressed genes
res <- results(dds, 
               contrast = c("condition", "PCOS", "control"))
## Order results by adjusted p-value
res <- res[order(res$padj), ]
head(res)
dim(res)

res_sig <- 
  res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 0.263),]
dim(res_sig)
tail(res_sig)
head(res_sig)


# b. Create MA plot
png("MAplot.png")
plotMA(res, main="DESeq2 MA Plot", ylim=c(-10,10))
abline(h=1, col="red", lwd=2, lty=2)
abline(h=-1, col="red", lwd=2, lty=2)
dev.off()


# c. # Volcano plot of the results
volcano_data <- as.data.frame(res)
volcano_data <- 
  volcano_data[!is.na(volcano_data$padj),]
volcano_data$significant <- 
  volcano_data$padj < 0.05 & abs(volcano_data$log2FoldChange) > 1

p <- ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(padj), 
                              color=significant)) +
  geom_point() +
  scale_color_manual(values=c("gray", "red")) +
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 Adjusted P-value") +
  theme_minimal()
print(p)
ggsave("volcano_plot.jpg", device = "jpeg")



# 9. Plot a heatmap of top 500 most significant genes
top_sig <- res_sig[1:500,]
top_sig <- rownames(top_sig)
head(top_sig)

# Get the matrix of top 500 genes with highest MADs
top_sig <- assay(vsd)[which(rownames(assay(vsd)) %in% top_sig),]
head(top_sig)

annot <- colData(dds)$condition
levels(annot)
category_colors <- c("control" = "blue", "PCOS" = "red")
mapped_colors <- category_colors[annot]

# Plot a heatmap of top 500 
png("top_signifcant_heatmap_clustering.png")
par(mar=c(2,2,2,2))
heatmap.2(top_sig,
          col=bluered(75),       
          trace="none",         
          dendrogram = "both",
          scale="row",          
          key=TRUE,              
          keysize=1.2,           
          margins=c(6,6),
          labRow = "",
          ColSideColors = mapped_colors
)

# Add a legend for the column annotations
legend("topright",
       legend=c("Normal", "PCOS"), 
       fill=category_colors, 
       title="Diagnosis",
       border=FALSE)
dev.off()





## Saving Results of DGE to use for GO and KEGG
# Convert DESeq2 results to a dataframe
res_df <- as.data.frame(res)  

# Preserve gene IDs
res_df$GeneID <- rownames(res_df) 

# Remove Gene_ID column
rownames(res_df) <- NULL  

# Saving the results in .csv
write.csv(res_df, file = "dge_results_unfiltered.csv", row.names = FALSE)

