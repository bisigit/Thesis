#############  Convert significant genes cpg sites to GeneID and Ensembl ID ##########
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)
library(ggpubr)


# Entrez
entrez_ids <- c("9973", "23250")


# 1. Convert to Ensembls
# Map to Ensembl IDs
ensembl_ids <- mapIds(
  org.Hs.eg.db,
  keys = entrez_ids,
  column = "ENSEMBL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Call results
ensembl_ids
# 9973             23250 
# "ENSG00000173992" "ENSG00000068650" 


# 2. Convert to cpg
# Step 1: Entrez → Gene Symbol
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = entrez_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# View symbols
head(gene_symbols)
# 9973    23250 
# "CCS" "ATP11A" 


# Step 2: Load annotation and find matching CpGs
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Filter for CpGs associated with those genes
matching_cpgs <- ann[grepl(paste(gene_symbols, collapse = "|"), ann$UCSC_RefGene_Name), ]

# View result
matching_cpgs[, c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]
cpgs <- matching_cpgs$Name


# Filter based on promoter cpgs
# Define promoter-related categories
promoter_contexts <- c("TSS200", "TSS1500", "5'UTR", "1stExon")

# Keep only rows where UCSC_RefGene_Group includes a promoter context
promoter_cpgs <- matching_cpgs[
  grepl(paste(promoter_contexts, collapse = "|"), matching_cpgs$UCSC_RefGene_Group),
]

# View the promoter CpGs and context
promoter_cpgs[, c("Name", "UCSC_RefGene_Name", "UCSC_RefGene_Group")]

# Extract just the CpG site names
cpgs_promoter_only <- promoter_cpgs$Name

print(cpgs_promoter_only)




########################## Getting boxplot for DMA ##########
# 1. Meth data
# Load
meth_matrix <- read.csv("beta_vals_for_DMA.csv", header = TRUE, row.names = 1)

# Subset to find cpgs that match 
cpgs_promoter_only <- intersect(row.names(meth_matrix), cpgs_promoter_only)

# Subset the beta matrix to include only those CpGs
beta_vals_subset <- meth_matrix[cpgs_promoter_only, , drop = FALSE]

# Check it worked
dim(beta_vals_subset)
head(beta_vals_subset)


# 2. RNA data
# Load
counts_matrix <- read.csv("cleaned_counts.csv", header = TRUE, row.names = 1)
metada <- read.csv("sample_info.csv", header = TRUE, row.names = 1)

# Create vector 
ensembl_ids_of_interest <- c("ENSG00000173992", "ENSG00000068650" )


# Filter by Ensembl IDs
# Intersect the Ensembl IDs with rownames
matched_ids <- intersect(rownames(counts_matrix), ensembl_ids_of_interest)

# Subset the expression matrix
filtered_counts <- counts_matrix[matched_ids, , drop = FALSE]

# View result
filtered_counts


# Running log transform counts
log_expr <- log2(filtered_counts + 1)



########################### Building boxplot using the above ###################
# 1. Build promoter table
# Make a simple table from promoter_cpgs
promoter_table_final <- data.frame(
  CpG = promoter_cpgs$Name,
  GeneSymbol = promoter_cpgs$UCSC_RefGene_Name,
  GeneContext = promoter_cpgs$UCSC_RefGene_Group,
  stringsAsFactors = FALSE
)

# Split rows where multiple genes are listed (e.g., "ATP11A;ATP11A")
library(tidyr)
promoter_table_final <- promoter_table_final %>%
  mutate(GeneSymbol = strsplit(GeneSymbol, ";")) %>%
  unnest(GeneSymbol)


# 2. Subset Cpgs per gene
# Get CpGs for each gene
atp11a_cpgs <- promoter_table_final %>% filter(GeneSymbol == "ATP11A") %>% pull(CpG)
ccs_cpgs    <- promoter_table_final %>% filter(GeneSymbol == "CCS") %>% pull(CpG)

# Compute mean beta per sample for each gene
atp11a_meth <- colMeans(beta_vals_subset[rownames(beta_vals_subset) %in% atp11a_cpgs, , drop = FALSE])
ccs_meth    <- colMeans(beta_vals_subset[rownames(beta_vals_subset) %in% ccs_cpgs,    , drop = FALSE])


# 3. Build dataframe for merging
# Match expression matrix rows to genes
ccs_expr    <- as.numeric(log_expr["ENSG00000173992", ])
atp11a_expr <- as.numeric(log_expr["ENSG00000068650", ])

# Final dataframe
cor_df <- data.frame(
  Sample = colnames(beta_vals_subset),
  CCS_Methylation = ccs_meth,
  ATP11A_Methylation = atp11a_meth,
  CCS_Expression = ccs_expr,
  ATP11A_Expression = atp11a_expr
)

# 4. Correlation scatter plot
# Reshape to long format
cor_df_long <- cor_df %>%
  select(Sample, CCS_Methylation, CCS_Expression, ATP11A_Methylation, ATP11A_Expression) %>%
  pivot_longer(cols = -Sample, names_to = c("Gene", ".value"), names_pattern = "(.*)_(Methylation|Expression)")

# Plot both genes in facets
ggplot(cor_df_long, aes(x = Methylation, y = Expression)) +
  geom_point(aes(color = Gene), size = 2, alpha = 0.8) +
  stat_cor(aes(color = Gene), method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4) +
  facet_wrap(~Gene, scales = "free") +
  labs(title = NULL,
       x = "Average Promoter Methylation", y = "Log Expression") +
  theme_minimal() +
  theme(legend.position = "none")
