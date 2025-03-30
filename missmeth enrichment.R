############### missmeth analysis ############

# Load libraries
library(missMethyl)
library(ggplot2)
library(ReactomePA)
library(org.Hs.eg.db)


# Set working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/METHarray data/Enrichment with missmeth")

# 1. Read the CSV you exported earlier
annotated_dmps <- read.csv("annotated_dmp_results.csv", row.names = 1)
betavals_for_dma <- read.csv("beta_vals_for_DMA.csv", row.names = 1)

# 2. Subset cpgs
# Significant CpGs
sig_cpgs <- annotated_dmps$TargetID

# All CpGs that were tested 
all_cpgs <- rownames(betavals_for_dma)


# 3. GO enrichment
go_results <- gometh(
  sig.cpg = sig_cpgs,
  all.cpg = all_cpgs,
  array.type = "450K",    
  collection = "GO"
)

# Vizualize
top_terms <- head(topGSA(go_results), 15)

go_plot <- ggplot(top_terms, aes(x = reorder(TERM, -P.DE), y = -log10(P.DE))) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = NULL, y = "-log10(p-value)", title = NULL) +
  theme_minimal()
ggsave("go_plot_latest.png", go_plot, width = 10, height = 6, dpi = 300)





# 4. KEGG enrichment
kegg_results <- gometh(
  sig.cpg = sig_cpgs,
  all.cpg = all_cpgs,
  array.type = "450K",
  collection = "KEGG"
)


# Vizualize
# Get top terms
KEGG_top <- head(topGSA(kegg_results), 15)

kegg_plot <- ggplot(KEGG_top, aes(x = reorder(Description, -P.DE), y = -log10(P.DE))) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = NULL, y = "-log10(p-value)", title = NULL) +
  theme_minimal()
ggsave("kegg_plot_latest.png", kegg_plot, width = 10, height = 6, dpi = 300)



# 5. Reactome enrichment
dmg_genes <- unique(unlist(strsplit(annotated_dmps$UCSC_RefGene_Name, ";|,")))

# Convert gene symbols to Entrez IDs (ReactomePA uses Entrez)
entrez_ids <- mapIds(org.Hs.eg.db, keys = dmg_genes,
                     column = "ENTREZID", keytype = "SYMBOL",
                     multiVals = "first")

# Run enrichment
reactome_res <- enrichPathway(gene = entrez_ids,
                              pvalueCutoff = 0.25,
                              readable = TRUE)

# View results
head(as.data.frame(reactome_res))

# Plot
png("reactome_barplot.png", width = 1200, height = 800)
barplot(reactome_res, showCategory = 15, title = "Reactome Pathway Enrichment")
dev.off()




# Saving the enrichment terms and lists from for analysis
