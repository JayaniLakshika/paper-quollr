library(Seurat)
library(Matrix)
library(phateR)
library(dplyr)
library(ggplot2)

# Set download path
download_path <- "~/Desktop/PhD Monash research files/Research papers/paper-quollr/data"
print(download_path)

# Download and extract data if it doesn't exist
if (!dir.exists(file.path(download_path, "scRNAseq", "T0_1A"))) {
  scprep::download_and_extract_zip(
    url = "https://md-datasets-public-files-prod.s3.eu-west-1.amazonaws.com/5739738f-d4dd-49f7-b8d1-5841abdbeb1e",
    output_dir = download_path
  )
}

# Load 10X data using Seurat v5
T1_seurat <- Read10X(data.dir = file.path(download_path, "scRNAseq", "T0_1A")) %>%
  CreateSeuratObject(project = "T0_1A")
T2_seurat <- Read10X(data.dir = file.path(download_path, "scRNAseq", "T2_3B")) %>%
  CreateSeuratObject(project = "T2_3B")
T3_seurat <- Read10X(data.dir = file.path(download_path, "scRNAseq", "T4_5C")) %>%
  CreateSeuratObject(project = "T4_5C")
T4_seurat <- Read10X(data.dir = file.path(download_path, "scRNAseq", "T6_7D")) %>%
  CreateSeuratObject(project = "T6_7D")
T5_seurat <- Read10X(data.dir = file.path(download_path, "scRNAseq", "T8_9E")) %>%
  CreateSeuratObject(project = "T8_9E")

# Combine Seurat objects
EBT_seurat <- merge(T1_seurat, y = c(T2_seurat, T3_seurat, T4_seurat, T5_seurat),
                    add.cell.ids = c("Day 00-03", "Day 06-09", "Day 12-15", "Day 18-21", "Day 24-27"))

# Filter cells based on library size (using Seurat's functions)
EBT_seurat <- subset(EBT_seurat, subset = nCount_RNA > quantile(EBT_seurat$nCount_RNA, 0.2))
EBT_seurat <- subset(EBT_seurat, subset = nCount_RNA < quantile(EBT_seurat$nCount_RNA, 0.75))

# Filter rare genes (using Seurat's functions)
counts_matrix <- GetAssayData(EBT_seurat, slot = "counts")
gene_counts <- rowSums(counts_matrix > 0)
genes_to_keep <- names(gene_counts[gene_counts >= 10])
EBT_seurat <- subset(EBT_seurat, features = genes_to_keep)

# Library size normalization (using Seurat's functions)
EBT_seurat <- NormalizeData(EBT_seurat)

# Get mitochondrial genes
mito_genes <- grep(pattern = "^MT-", x = rownames(EBT_se-urat), value = TRUE)

# Filter cells based on mitochondrial gene expression
EBT_seurat <- subset(EBT_seurat, subset = percent.mt < quantile(EBT_seurat$percent.mt, 0.9))

# Square root transformation (using Seurat's functions)
EBT_seurat <- NormalizeData(EBT_seurat, normalization.method = "CLR")
EBT_seurat[["RNA"]] <- ScaleData(EBT_seurat, features = rownames(EBT_seurat), normalization.method = "center")

# Access the count matrix
EBT_counts <- GetAssayData(EBT_seurat, slot = "counts")
sample_labels <- EBT_seurat$orig.ident

# Display the head of EBT_counts
print(head(EBT_counts))
