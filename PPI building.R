########## Building a PPI Network using the genes from top pathways enriched 

# Load libraries
library(STRINGdb)
library(igraph)
library(org.Hs.eg.db)
library(biomaRt)



# Set working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/RNAseq data/Reactome Go and KEGG/PPI Network")


# 1. Load data
genes_enriched <- read.csv("ppi_leading_genes_with_expression.csv")

# Confirm EntrezID:
head(genes_enriched$EntrezID)


# 2. Create object
options(timeout = 300)
string_db <- STRINGdb$new(version = "11.5",
                          species = 9606,
                          score_threshold = 700,  # high confidence
                          input_directory = "")

# 3. Map genes to String IDs
# Add genes as a data frame with column "entrez"
genes_enriched$gene <- as.character(genes_enriched$EntrezID)

# Map to STRING
mapped_genes <- string_db$map(genes_enriched,
                              my_data_frame_id_col_names = "gene",  
                              removeUnmappedRows = TRUE)

# 4. Get all high-confidence interactions
interactions <- string_db$get_interactions(mapped_genes$STRING_id)





# Get Entrez → Symbol mapping
gene_symbols <- bitr(mapped_genes$gene,
                     fromType = "ENTREZID",
                     toType = "SYMBOL",
                     OrgDb = org.Hs.eg.db)

# 5. Add SYMBOL to mapped_genes for merging
mapped_genes <- merge(mapped_genes, gene_symbols,
                      by.x = "gene", by.y = "ENTREZID")

# 6. Prepare node labels
node_labels <- mapped_genes[, c("STRING_id", "SYMBOL")]

# 7. Merge with interactions (from and to)
interactions_named <- merge(interactions, node_labels,
                            by.x = "from", by.y = "STRING_id", all.x = TRUE)
colnames(interactions_named)[ncol(interactions_named)] <- "from_symbol"

interactions_named <- merge(interactions_named, node_labels,
                            by.x = "to", by.y = "STRING_id", all.x = TRUE)
colnames(interactions_named)[ncol(interactions_named)] <- "to_symbol"

# 8.Keep only clean columns for Cytoscape
ppi_for_cytoscape <- interactions_named[, c("from_symbol", "to_symbol", "combined_score")]

# 9. Save
write.csv(ppi_for_cytoscape, "ppi_for_cytoscape.csv", row.names = FALSE, quote = FALSE)


# 10. Create node table for Cytoscape
node_table <- mapped_genes[, c("SYMBOL", "log2FoldChange", "pvalue", "padj")]

# Remove duplicate rows if any
node_table <- node_table[!duplicated(node_table$SYMBOL), ]

# Remove quotes from SYMBOLs
node_table$SYMBOL <- gsub('"', '', node_table$SYMBOL)

# Just to be safe: also trim whitespace
node_table$SYMBOL <- trimws(node_table$SYMBOL)

# Save again
write.csv(node_table, "cytoscape_node_attributes_clean.csv", row.names = FALSE, quote = FALSE)

