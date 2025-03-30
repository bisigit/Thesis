############## Overlap of DGE and DMRS ######################

# Load libraries
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(org.Hs.eg.db)
library(ggVennDiagram)
library(ggplot2)


# Set Working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/METHarray data/Enrichment with missmeth")


# 1. Load DEGs (with Entrez IDs)
deg_df <- read.csv("dge_sig_res_entrezID.csv")   
deg_entrez <- unique(deg_df$entrez)  

# Load DMRs (with CpG TargetIDs)
dmr_df <- read.csv("dmp_results latest.csv")  
dmr_cpgs <- unique(dmr_df$TargetID)    

# Load CTD gene list
ctd_df <- read.csv("ctd_pcos_genes.csv")  
ctd_entrez <- unique(ctd_df$Gene.ID) 

# 2. Annotate CpG IDs to Entrez IDs
# Load annotation for 450K array
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Map CpG probe IDs to gene symbols
dmr_annotated <- merge(dmr_df, ann, by.x = "TargetID", by.y = "Name")

# Extract gene symbols from CpGs
dmr_genes <- unique(unlist(strsplit(dmr_annotated$UCSC_RefGene_Name, ";|,")))

# Map gene symbols to Entrez IDs
dmr_entrez <- mapIds(org.Hs.eg.db, keys = dmr_genes,
                     column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
dmr_entrez <- na.omit(unique(dmr_entrez))
write.csv(dmr_entrez, "dmr_entrez.csv")


# 3. Find overlap
overlap_genes <- intersect(deg_entrez, dmr_entrez)

# Report
cat("Number of overlapping genes:", length(overlap_genes), "\n")
print(overlap_genes)
# [1] "9973"  "23250"



# 4. Find out what the overlapping genes are
# Map Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = overlap_genes,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")


# View results
print(gene_symbols)
# 9973    23250 
# "CCS" "ATP11A" 



########################### Overlap of DEGs and DMRs with known PCOS genes ##########################
# 1. Overlap with DEGs
# Run overlap
deg_ctd_overlap <- intersect(deg_entrez, ctd_entrez)  # 2,350 overlaps
write.csv(deg_ctd_overlap, "deg_ctd_overlap_entrez.csv")

# DEG overlap gene symbols
deg_symbols <- mapIds(org.Hs.eg.db, keys = as.character(deg_ctd_overlap),
                      column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
write.csv(deg_symbols, "deg_ctd_overlap_symbols.csv")


# Overlap with DMRs
dmr_ctd_overlap <- intersect(dmr_entrez, ctd_entrez)  # 17 overlaps
write.csv(dmr_ctd_overlap, "dmr_ctd_overlap_entrez.csv")

# DMR overlap gene symbols
dmr_symbols <- mapIds(org.Hs.eg.db, keys = dmr_ctd_overlap,
                      column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
write.csv(dmr_symbols, "dmr_ctd_overlap_symbols.csv")




########################### Visualize Overlaps with Venn Diagrams ################
# 1. DMR, DGE, CTD Venn Diagram
# Create named list
gene_sets <- list(
  DEGs = as.character(deg_entrez),
  DMRs = as.character(dmr_entrez),
  CTD_PCOS = as.character(ctd_entrez)
)

# Draw the Venn Diagram
vennd <- ggVennDiagram(gene_sets, label = "count") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal()
ggsave("dge_dmr_ctd_venn.png", plot = vennd)



# 2. DMR and DGE Venn Diagram
# Create a new named list with just DEGs and DMRs
gene_sets_2way <- list(
  DEGs = as.character(deg_entrez),
  DMRs = as.character(dmr_entrez)
)

# Draw the 2-way Venn diagram
venn_2way <- ggVennDiagram(gene_sets_2way, label = "count") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  # Fix background
    plot.background = element_rect(fill = "white", color = NA)    # Fix plot background
  )
ggsave("dge_dmr_venn.png", plot = venn_2way, width = 6, height = 6, dpi = 300)
