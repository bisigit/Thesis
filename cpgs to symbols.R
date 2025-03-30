############## converting methylation cpgs to gene symbols for PPI#########

# Load libraries
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(STRINGdb)
library(igraph)



# Set working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/METHarray data/PPI building")


# 1. Load data
meth_results <- read.csv("dmp_results latest.csv", row.names = 1)

# 2. Map cpgs
# Load annotation
annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Extract matching rows from the annotation
annot_subset <- annotation[rownames(meth_results), ]

# Add gene symbol to methylation results
meth_results$Gene <- annot_subset$UCSC_RefGene_Name

# Split results with many mappings
meth_results_expanded <- meth_results %>%
  separate_rows(Gene, sep = ";|,")







################## Mapping to stringdb ##################
# 1. Deduplicate gene column
unique_genes_df <- meth_results_expanded %>%
  filter(!is.na(Gene), Gene != "") %>%  
  distinct(Gene, .keep_all = FALSE)


# Remove blanks
genes_dmg <- data.frame(gene = unique_genes_df$Gene)


# 2. Map to Stringdb
# Create STRINGdb object
string_db <- STRINGdb$new(version = "11.5",
                          species = 9606,
                          score_threshold = 200,
                          input_directory = "")

# Map gene symbols to STRING IDs
mapped_genes_dmg <- string_db$map(genes_dmg,
                                  my_data_frame_id_col_names = "gene",
                                  removeUnmappedRows = TRUE)

valid_ids <- unique(as.character(na.omit(mapped_genes_dmg$STRING_id)))

# Check it manually
length(valid_ids)  # how many unique, non-NA IDs
head(valid_ids)

# Retrieve interactions
options(timeout = 300)  
interactions_dmg <- string_db$get_interactions(valid_ids)


# 3. Map back to gene symbols 
# Merge back with mapped gene symbols
node_labels <- mapped_genes_dmg[, c("STRING_id", "gene")]

# Join to get from_symbol and to_symbol
interactions_named <- merge(interactions_dmg, node_labels,
                            by.x = "from", by.y = "STRING_id", all.x = TRUE)
colnames(interactions_named)[ncol(interactions_named)] <- "from_symbol"

interactions_named <- merge(interactions_named, node_labels,
                            by.x = "to", by.y = "STRING_id", all.x = TRUE)
colnames(interactions_named)[ncol(interactions_named)] <- "to_symbol"

# Optional: keep only complete rows
ppi_for_cytoscape_dmg <- interactions_named[, c("from_symbol", "to_symbol", "combined_score")]
ppi_for_cytoscape_dmg <- na.omit(ppi_for_cytoscape_dmg)

# Save file for cytoscape
write.csv(ppi_for_cytoscape_dmg, "dmr_data_for_cytoscape.csv", row.names = FALSE, quote = FALSE)


# 4. Get file with logfc values for expression data in cytoscape
# Get genes involved in PPI
genes_in_ppi <- unique(c(ppi_for_cytoscape_dmg$from_symbol,
                         ppi_for_cytoscape_dmg$to_symbol))

# Remove NAs or blanks just in case
genes_in_ppi <- genes_in_ppi[!is.na(genes_in_ppi) & genes_in_ppi != ""]

# Subset methylation results
meth_for_cytoscape <- meth_results_expanded %>%
  filter(Gene %in% genes_in_ppi)

# Keep only columns for Cytoscape styling
node_table <- meth_for_cytoscape %>%
  select(Gene, logFC, adj.P.Val) %>%
  distinct(Gene, .keep_all = TRUE)  

# For cytoscape
write.csv(node_table, "dmg_node_attributes.csv", row.names = FALSE, quote = FALSE)

