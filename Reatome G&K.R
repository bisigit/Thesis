################ Running GSEA with Reactome ##################

# Load Libraries
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(patchwork)

# Set working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/RNAseq data/Reactome Go and KEGG")


# 1. Set Up: Load data
# Import a CSV file into R
results <- fread("dge_sig_results_filtered.csv")

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
write.csv(results, "dge_sig_res_entrezID.csv", row.names = TRUE)




# 2. GSEA with Reactome: Pathway analysis 
# gene_ranks: named vector of log2FC values with Entrez IDs as names
gene_ranks <- results$log2FoldChange
names(gene_ranks) <- results$entrez  

# Remove NAs
gene_ranks <- gene_ranks[!is.na(names(gene_ranks)) & !is.na(gene_ranks)]

# Sort in decreasing order.
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# Run GSEA
reactome_gsea <- gsePathway(geneList = gene_ranks,
                            organism = "human",
                            pvalueCutoff = 0.25,
                            verbose = FALSE)

# Plot Top 15 Reactome Pathways - dotplot
p_reactome <- dotplot(reactome_gsea, showCategory = 15, title = NULL)
ggsave("kegg_reactome_dotplot.pdf", plot = p_reactome,
       width = 8, height = 6, dpi = 300)




# 3. GSEA for GO:BP, GO:CC, and GO:MF
# GO:BP (Biological Process)
gsea_bp <- gseGO(geneList = gene_ranks,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 keyType = "ENTREZID",
                 pvalueCutoff = 0.25)

# GO:CC (Cellular Component)
gsea_cc <- gseGO(geneList = gene_ranks,
                 OrgDb = org.Hs.eg.db,
                 ont = "CC",
                 keyType = "ENTREZID",
                 pvalueCutoff = 0.25)

# GO:MF (Molecular Function)
gsea_mf <- gseGO(geneList = gene_ranks,
                 OrgDb = org.Hs.eg.db,
                 ont = "MF",
                 keyType = "ENTREZID",
                 pvalueCutoff = 0.25)


# Vizualize
# Biological Process
p_bp <- dotplot(gsea_bp, showCategory = 15, title = "GSEA GO: Biological Process") +
  theme(axis.text.y = element_text(size = 12),      
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))

# Cellular Component
p_cc <- dotplot(gsea_cc, showCategory = 15, title = "GSEA GO: Cellular Component") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))

# Molecular Function
p_mf <- dotplot(gsea_mf, showCategory = 15, title = "GSEA GO: Molecular Function") +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))

# Save plots
ggsave("gsea_go_bp_dotplot.pdf", plot = p_bp, width = 10, height = 6, dpi = 300)
ggsave("gsea_go_cc_dotplot.pdf", plot = p_cc, width = 10, height = 6, dpi = 300)
ggsave("gsea_go_mf_dotplot.pdf", plot = p_mf, width = 10, height = 6, dpi = 300)





# Additional Data to export
# 1. Save results to CSV
write.csv(gsea_bp@result, "GSEA_GO_BP_full_results.csv", row.names = FALSE)
write.csv(gsea_cc@result, "GSEA_GO_CC_full_results.csv", row.names = FALSE)
write.csv(gsea_mf@result, "GSEA_GO_MF_full_results.csv", row.names = FALSE)
write.csv(reactome_gsea@result, "GSEA_Reactome_full_results.csv", row.names = FALSE)



# 2. Get leading-edge genes for
# Specify desired GO terms
target_bp <- c("protein phosphorylation", "immune response", "carbohydrate metabolic process")
target_cc <- c("membrane raft", "plasma membrane protein complex")
target_mf <- c("protein kinase activity", "signaling receptor binding")
target_reactome <- c("TP53 Regulates Transcription of DNA Repair Genes", "Immune System / Innate Immune System", "Interferon alpha/beta signaling")

# Extract all leading-edge genes from selected GO terms
extract_all_leading_genes <- function(gsea_result, target_terms) {
  genes <- gsea_result@result %>%
    dplyr::filter(Description %in% target_terms) %>%
    dplyr::pull(core_enrichment) %>%
    strsplit("/") %>%
    unlist()
  return(genes)
}

# Extract genes
genes_bp <- extract_all_leading_genes(gsea_bp, target_bp)
genes_cc <- extract_all_leading_genes(gsea_cc, target_cc)
genes_mf <- extract_all_leading_genes(gsea_mf, target_mf)
genes_reactome <- extract_all_leading_genes(reactome_gsea, target_reactome)

# Combine genes
all_leading_genes <- c(genes_bp, genes_cc, genes_mf, genes_reactome)

# Count frequency of each gene's appearance across terms
library(dplyr)
gene_counts <- as.data.frame(table(all_leading_genes))
gene_counts <- gene_counts %>%
  arrange(desc(Freq)) %>%
  slice_head(n = 100)  


# Combine all
ppi_gene_list <- gene_counts$all_leading_genes

# Convert to data frame for export
leading_genes_df <- data.frame(EntrezID = ppi_gene_list)

# Merge leading-edge genes with full expression data
leading_expr_data <- merge(leading_genes_df, results, by.x = "EntrezID", by.y = "entrez", all.x = TRUE)

# Save merged data to CSV
write.csv(leading_expr_data, "ppi_leading_genes_with_expression.csv", row.names = FALSE)







# 3. Category network plot --- I'M NOT USING THIS ANYMORE
# a. Biological Process
bp_net <- cnetplot(gsea_bp, showCategory = 5, foldChange = gene_ranks) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
ggsave("networkplot_bp.png", plot = bp_net, width = 10, height = 6, dpi = 300)

# b. Cellular Component
cc_net <- cnetplot(gsea_cc, showCategory = 5, foldChange = gene_ranks) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
ggsave("networkplot_cc.png", plot = cc_net, width = 10, height = 6, dpi = 300)

# c. Molecular Function
mf_net <- cnetplot(gsea_mf, showCategory = 5, foldChange = gene_ranks) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))
ggsave("networkplot_mf.png", plot = mf_net, width = 10, height = 6, dpi = 300)

# Combined plot
# If not installed yet
# install.packages("patchwork")

library(patchwork)

# Combine vertically (you can also use `|` for horizontal)
combined_net <- bp_net / cc_net / mf_net +
  plot_annotation(title = "GO Category Network Plots",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

# Save
ggsave("networkplot_combined.png", plot = combined_net,
       width = 12, height = 18, dpi = 300)
