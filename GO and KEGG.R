## GO and KEGG on GSE84958 ##

# Loading relevant libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(biomaRt)
library(ggplot2)
library(readr)
library(data.table)

# Set working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/RNAseq data")


# Load data
# Import a CSV file into R
results <- fread("dge_results_unfiltered.csv")

# Set row names as GeneID
rownames(results) <- results$GeneID

# Delete GeneID column
results$GeneID <- NULL

# Mapping Entrez IDs for results (res)
results$entrez <- mapIds(org.Hs.eg.db, 
                        keys = rownames(results), 
                        column = "ENTREZID", 
                        keytype = "ENSEMBL", 
                        multiVals = "first")

# Removing rows where Entrez ID is missing 
results <- results[!is.na(results$entrez), ]

# Remove duplicate Entrez IDs by keeping the highest absolute log2FoldChange
results <- results[!duplicated(results$entrez), ] 


# I. GO
# 1. Run Over Representation Analysis (enrichGO)
# Filter significant DEGs 
results_sig <- 
  results[which(results$padj < 0.05 & abs(results$log2FoldChange) > 0.263),]

# Run ORA
ego <- enrichGO(
  gene = results_sig$entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE  
)

# Dotplot
p1 <- dotplot(ego, showCategory=10)
ggsave("ora_dotplot.png", plot = p1, width = 8, height = 6, dpi = 300)

# Barplot
p2 <- barplot(ego, showCategory=10)
ggsave("ora_barplot.png", plot = p2, width = 8, height = 6, dpi = 300)


# 2. Run gene set enrichment analysis (gseGEO)
gene_ranks <- results_sig$log2FoldChange
names(gene_ranks) <- results_sig$entrez
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# Run gseGO
gsea_res <- gseGO(
  geneList = gene_ranks,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP", 
  pvalueCutoff = 0.05
)

# Ridgeplot
p11 <- ridgeplot(gsea_res)
ggsave("GSEA_ridgeplot.png", plot = p11, width = 8, height = 6, dpi = 300)

# Enrichment Plot
p12 <- gseaplot2(gsea_res, geneSetID=1, title="Top GO Term")
ggsave("GSEA_plot2.png", plot = p12, width = 8, height = 6, dpi = 300)

# Dotplot
pic1 <- dotplot(ego, showCategory=10)
ggsave("GSEA_dotplot.png", plot = pic1, width = 8, height = 6, dpi = 300)

# Barplot
pic2 <- barplot(ego, showCategory=10)
ggsave("GSEA_barplot.png", plot = pic2, width = 8, height = 6, dpi = 300)



# II. KEGG
# 1. Over Representation Analysis
# Run enrichKEGG
kegg_ora <- enrichKEGG(
  gene = results_sig$entrez,
  organism = "hsa",  
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# Convert Entrez IDs back to gene symbols for readability
kegg_ora <- setReadable(kegg_ora, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Dotplot
dotplot(kegg_ora, showCategory=10)
ggsave("KEGG_ora_dotplot.png", width = 8, height = 6, dpi = 300)

# Barplot
barplot(kegg_ora, showCategory=10)
ggsave("KEGG_ora_barplot.png", width = 8, height = 6, dpi = 300)


# 2. Gene Set Enrichment analysis for KEGG
kegg_gsea <- gseKEGG(
  geneList = gene_ranks,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# Ridgeplot
ridgeplot(kegg_gsea)
ggsave("KEGG_gsea_ridgeplot.png", width = 8, height = 6, dpi = 300)

# GSEA Enrichment plot for top KEGG pathway
gseaplot2(kegg_gsea, geneSetID=1, title="Top KEGG Pathway")
ggsave("KEGG_gsea_top_pathway.png", width = 8, height = 6, dpi = 300)


# Dotplot
pic11 <- dotplot(kegg_gsea, showCategory=10)
ggsave("KEGG_gsea_dotplot.png", plot = pic11, width = 8, height = 6, dpi = 300)


