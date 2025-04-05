## Tabula Muris: single organ analysis (Limb Muscles)
## https://github.com/czbiohub-sf/tabula-muris/blob/master/00_data_ingest/FACS_Notebook.Rmd

library(Seurat)
library(dplyr)
library(Matrix)
library(stringr)
library(here)
library(tidyverse)

organ = "Limb_Muscle"

raw.data <- read.csv(here("data", "FACS",paste0(organ,"-counts.csv")),row.names = 1)
meta.data <- read.csv(here("data", "metadata_FACS.csv"))

plates <- str_split(colnames(raw.data),"[.]", simplify = TRUE)[,2]

rownames(meta.data) <- meta.data$plate.barcode
cell.meta.data <- meta.data[plates,]
rownames(cell.meta.data) <- colnames(raw.data)

# Find ERCC's, compute the percent ERCC, and drop them from the raw data.
erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
raw.data <- raw.data[-ercc.index,]

# Create the Seurat object with all the data
tiss <- CreateSeuratObject(counts = raw.data, min.cells = 5, min.features = 5) # Changed raw.data to counts

tiss <- AddMetaData(object = tiss, cell.meta.data)
tiss[["percent.ercc"]] <- percent.ercc # Updated AddMetaData for v5
# Change default name for sums of counts from nUMI to nReads
tiss@meta.data <- tiss@meta.data %>% rename(nReads = nCount_RNA, nGene = nFeature_RNA) # Updated rename for v5

ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = tiss), value = TRUE) # Changed tiss@data to tiss
percent.ribo <- Matrix::colSums(GetAssayData(tiss, slot = "counts")[ribo.genes, ])/Matrix::colSums(GetAssayData(tiss, slot = "counts"))
tiss[["percent.ribo"]] <- percent.ribo # Updated AddMetaData for v5

percent.Rn45s <- Matrix::colSums(as.matrix(GetAssayData(tiss, slot = "counts")['Rn45s', , drop = FALSE])) / Matrix::colSums(GetAssayData(tiss, slot = "counts"))
tiss[["percent.Rn45s"]] <- percent.Rn45s # Updated AddMetaData for v5

FeatureScatter(object = tiss, feature1 = "nReads", feature2 = "nGene") # Changed GenePlot to FeatureScatter and use [email address removed] columns

tiss <- subset(tiss, subset = nGene > 500 & nReads > 50000 & nGene < 25000 & nReads < 2000000) # Updated FilterCells to subset

tiss <- NormalizeData(object = tiss)
tiss <- ScaleData(object = tiss, vars.to.regress = c("nReads", "percent.ribo","percent.Rn45s")) # Updated vars.to.regress to match metadata
tiss <- FindVariableFeatures(object = tiss, selection.method = "vst", nfeatures = 2000) # Updated FindVariableGenes to FindVariableFeatures

tiss <- RunPCA(object = tiss, verbose = FALSE) # Updated do.print to verbose
tiss <- ProjectDim(object = tiss, verbose = FALSE) # Updated ProjectPCA to ProjectDim

ElbowPlot(object = tiss, ndims = 3)

pca_df <- tiss@reductions$pca@cell.embeddings
pca_df <- pca_df[, 1:10] |>
  tibble::as_tibble()

names(pca_df) <- paste0("x", 1:10)

write_rds(pca_df, "data/limb_muscles/facs_limb_muscles_pcs_10.rds")
