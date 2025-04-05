## Tabula Muris: single organ analysis (Limb Muscles)
## https://github.com/czbiohub-sf/tabula-muris/blob/master/00_data_ingest/FACS_Notebook.Rmd

library(Seurat)
library(dplyr)
library(Matrix)
library(stringr)
library(here)

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
tiss <- CreateSeuratObject(raw.data = raw.data, min.cells = 5, min.genes = 5)

tiss <- AddMetaData(object = tiss, cell.meta.data)
tiss <- AddMetaData(object = tiss, percent.ercc, col.name = "percent.ercc")
# Change default name for sums of counts from nUMI to nReads
colnames(tiss@meta.data)[colnames(tiss@meta.data) == 'nUMI'] <- 'nReads'

ribo.genes <- grep(pattern = "^Rp[sl][[:digit:]]", x = rownames(x = tiss@data), value = TRUE)
percent.ribo <- Matrix::colSums(tiss@raw.data[ribo.genes, ])/Matrix::colSums(tiss@raw.data)
tiss <- AddMetaData(object = tiss, metadata = percent.ribo, col.name = "percent.ribo")

percent.Rn45s <- Matrix::colSums(tiss@raw.data[c('Rn45s'), ])/Matrix::colSums(tiss@raw.data)
tiss <- AddMetaData(object = tiss, metadata = percent.Rn45s, col.name = "percent.Rn45s")

GenePlot(object = tiss, gene1 = "nReads", gene2 = "nGene", use.raw=T)

tiss <- FilterCells(object = tiss, subset.names = c("nGene", "nReads"),
                    low.thresholds = c(500, 50000), high.thresholds = c(25000, 2000000))

tiss <- NormalizeData(object = tiss)
tiss <- ScaleData(object = tiss, vars.to.regress = c("nReads", "percent.ribo","Rn45s"))
tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5)

tiss <- RunPCA(object = tiss, do.print = FALSE)
tiss <- ProjectPCA(object = tiss, do.print = FALSE)

PCHeatmap(object = tiss, pc.use = 1:3, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, num.genes = 8)

PCElbowPlot(object = tiss)

# Set number of principal components.
n.pcs = 10


