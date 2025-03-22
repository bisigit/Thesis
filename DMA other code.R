## DMA with Beta Values ##

# Load libraries
library(readr)
library(data.table)
library(dplyr)
library(tidyverse)
library(minfi)
library(limma)
library(impute)
library(EnhancedVolcano)
library(GEOquery)
library(sva)
library(ggplot2)
library(preprocessCore)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Set working directory
setwd("C:/Users/amanw/OneDrive - University of Lethbridge/Desktop/Winter 2025/Thesis/METHarray data")

# 1. Manually download dataset file the files section in ArrayExpress

# 2. Read the processed methylation data
meth_a <- fread("PCOS_SubQ_Adipose_Methylation_450k_a.txt", header = TRUE)
meth_b <- fread("PCOS_SubQ_Adipose_Methylation_450k_b.txt", header = TRUE)
metadata <- fread("E-MTAB-3777.sdrf.txt", header = TRUE)

# Convert to data frame
meth_a <- as.data.frame(meth_a)
meth_b <- as.data.frame(meth_b)

# Check if CpG IDs match between files
all(rownames(meth_a) == rownames(meth_b))  

# Check for duplicate sample names
colnames(meth_b)[!(colnames(meth_b) %in% colnames(meth_a))]  ##"ILMNID", "NAME", and "V114"  looked pretty different from the other columns

# Check if "ILMNID", "NAME", and "V114" contain redundant information
head(meth_b[, c("ILMNID", "NAME", "V114")]) #Yes

# Check the differences in column names between the two samples
setdiff(colnames(meth_a), colnames(meth_b))   # column v112 is unique
setdiff(colnames(meth_b), colnames(meth_a))   # ILMNID and Names are redundant but column V114 is unique

# Check content of V112 and V114
meth_a$V112  # contains N/A values
meth_b$V114  # contains N/A

# Check if all the values in v112 and v114 are N/A
all(is.na(meth_a$V112))  # TRUE
all(is.na(meth_b$V114))  # TRUE

# Remove column V112 from meth_a
meth_a <- meth_a %>% select(-V112)

# Remove the 3 redundant columns from meth_b
meth_b <- meth_b %>% select(-c(ILMNID, NAME, V114))

# 3. Merge the dataframes by TargetID, probeID_A, and probeID_B
merged_data <- inner_join(meth_a, meth_b, by = c("TargetID", "ProbeID_A", "ProbeID_B"))

# Set the TargetID (CpG sites) as row names 
rownames(merged_data) <- meth_a$TargetID

# Remove the TargetID 
merged_data$TargetID <- NULL

# 4. Getting the CpG x sample matrix
# Extract sample IDs
column_names <- colnames(merged_data)
sample_ids <- unique(str_extract(column_names, "^\\d+"))
sample_ids <- sample_ids[!is.na(sample_ids)]

# Extract AVG_Beta Columns and Tidy
tidy_beta <- merged_data %>%
  mutate(cpg_site = rownames(merged_data)) %>%
  pivot_longer(cols = -cpg_site, names_to = "variable", values_to = "AVG_Beta") %>%
  filter(str_detect(variable, "AVG_Beta")) %>%
  mutate(sample_id = str_extract(variable, "^\\d+")) %>%
  select(cpg_site, sample_id, AVG_Beta)

# Create a wide dataframe.
wide_beta <- tidy_beta %>%
  pivot_wider(names_from = sample_id, values_from = AVG_Beta)

# Convert to dataframe
wide_beta <- as.data.frame(wide_beta)  

# Reset the TargetID (CpG sites) as row names 
rownames(wide_beta) <- wide_beta$cpg_site

# Remove the TargetID 
wide_beta$cpg_site <- NULL

# 5. Save sample-by-CpG matrix for analyses
write.csv(wide_beta, "meth_matrix.csv", row.names = TRUE)

wide_beta <- fread("meth_matrix.csv", header = TRUE)
rownames(wide_beta) <- wide_beta$V1
wide_beta$V1 <- NULL
metadata <- fread("E-MTAB-3777.sdrf.txt", header = TRUE)

######################### Processing Data for DMA ##############################

# 1. Normalize beta values
# Check before normalization
boxplot(wide_beta, main = "B-values Before Normalization", outline = FALSE)

# Normalization
beta_values <- normalize.quantiles(as.matrix(wide_beta))

# Check after normalization
boxplot(beta_values, main = "B-values After Normalization", outline = FALSE)

# Reset column and row names
colnames(beta_values) <- colnames(wide_beta)
row.names(beta_values) <- row.names(wide_beta)

# Check for NA values
sum(is.na(beta_values))

# Calculate the percentage of NA values
na_percentage <- sum(is.na(beta_values)) / (nrow(beta_values) * ncol(beta_values)) * 100
print(paste("Percentage of NA values:", na_percentage, "%"))    # [1] "Percentage of NA values: 0.14476763383219 %"

# 2. Remove rows with fix missing values
complete_data <- beta_values[complete.cases(beta_values), ]

# 3. Extract sample names from metadata and changing columnnames to disease state
sample_names <- metadata[, c("Source Name", "Characteristics[disease]")]
colnames(sample_names) <- c("SampleName", "DiseaseState")

# Convert disease state names
sample_names$DiseaseState <- recode(sample_names$DiseaseState, 
                                    "normal" = "control", 
                                    "polycystic ovary syndrome" = "PCOS")

# Match sample order in methylation matrix with metadata
sample_names <- sample_names[match(colnames(complete_data), sample_names$SampleName), , drop = FALSE]

# Rename the columns of normalized_m_values based on 'Characteristics[disease]'
colnames(complete_data) <- make.unique(sample_names$`DiseaseState`)

# Ensure complete_data is a data frame
complete_data <- as.data.frame(complete_data)

# 4. # Remove probes with extreme beta values (all 0s or all 1s)
keep_probes <- rowSums(complete_data == 0) < (ncol(complete_data) * 0.95) &
  rowSums(complete_data == 1) < (ncol(complete_data) * 0.95)
complete_data <- complete_data[keep_probes, ]

# 4. PCA 
# Run PCA
pca_result <- prcomp(t(complete_data), center = TRUE, scale. = TRUE)

# PCA plot (scatter plot of the first two principal components)
pca_data <- data.frame(pca_result$x)
pca_data$SampleID <- colnames(complete_data)

# Scatter plot of the first two principal components
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, label = SampleID)) +
  geom_point(size = 3) +
  geom_text(aes(label = SampleID), vjust = -1, size = 3) +
  labs(title = "PCA of Methylation Data", x = "PC1", y = "PC2") +
  theme_minimal()
ggsave("pca_plot.png", plot = pca_plot, width = 8, height = 6)  #outlier: PCOS.2

# 5. Remove the outlier column
#complete_data$`PCOS.2` <- NULL

# Remove the row in sample_names corresponding to the outlier
#sample_names <- sample_names[-5, ]
#sample_names <- sample_names[-which(sample_names$SampleName=="PCOS.2"),]

# 6. SVA analysis
# Defining model matrices
mod <- model.matrix(~ 0 + DiseaseState, data = sample_names)
mod0 <- model.matrix(~1, sample_names)

# Run sva
# Estimate the number of surrogate variables
num.sv <- num.sv(complete_data, mod, method = "be") # 6

# Run SVA to obtain surrogate variables
sva_results <- sva(as.matrix(complete_data), mod, mod0, n.sv = num.sv)

# Get surrogate variables
surrogate <- sva_results$sv
surrogate <- surrogate[,1:6]
colnames(surrogate) <- c("SV1", "SV2", "SV3", "SV4", "SV5",
                         "SV6")

# Add surrogate variables to sample data
sample_names <- cbind(sample_names, surrogate) 

# # Ensure disease state is a factor with levels
sample_names$DiseaseState <- factor(sample_names$DiseaseState, levels = c("control", "PCOS"))
# Create model matrix with adjustment and surrogate variables
design <- model.matrix(
  ~ DiseaseState + SV1 +
    SV2 + SV3 + SV4 + SV5 + SV6,
  data = sample_names
)

# 7. DMR analysis with Limma
# Convert methylation data to a matrix (important for limma)
completed_data <- as.matrix(complete_data)

# Fit linear model
fit <- lmFit(completed_data, design)

# Apply contrasts and empirical Bayes moderation
fit2 <- eBayes(fit)

# Extract differentially methylated probes (DMPs)
dmps <- topTable(fit2, coef = "DiseaseStatePCOS", number = Inf, adjust.method = "BH", sort.by = "P")

