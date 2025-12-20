library(dplyr)
library(tibble)
library(readr)
library(conflicted)

library(Rtsne)
library(umap)
library(phateR)
library(reticulate)

set.seed(20240110)

conflicts_prefer(umap::umap)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

use_python("~/miniforge3/envs/pcamp_env/bin/python")
use_condaenv("pcamp_env")

data <- read_rds(here::here("data/limb_muscles/facs_limb_muscles_pcs_10.rds"))
data <- data |>
  dplyr::select(-ID)

## tSNE (default)
perplexity <- 30

tSNE_fit <- data |>
  dplyr::select(where(is.numeric)) |>
  Rtsne::Rtsne(perplexity = perplexity,
               pca = FALSE)

tSNE_data <- tSNE_fit$Y |>
  tibble::as_tibble(.name_repair = "unique")
names(tSNE_data) <- c("emb1", "emb2")

tSNE_data <- tSNE_data |>
  mutate(ID = row_number())

write_rds(tSNE_data, file = paste0("data/limb_muscles/facs_limb_muscles_tsne_perplexity_", perplexity, ".rds"))

## UMAP

# Create a config list with the desired parameters
umap_config <- umap.defaults
umap_config$n_neighbors <- 15      # Set the number of neighbors
umap_config$n_components <- 2    # Set the number of output dimensions (typically 2 or 3)
umap_config$min_dist <- 0.1

UMAP_fit <- umap(data, config = umap_config)

UMAP_data <- UMAP_fit$layout |>
  as_tibble()

names(UMAP_data) <- c("emb1", "emb2")

UMAP_data <- UMAP_data |>
  mutate(ID = row_number())

## Run only once
write_rds(UMAP_data, file = paste0("data/limb_muscles/facs_limb_muscles_umap_n-neigbors_", "15", "_min-dist_", "0.1", ".rds"))

## PHATE
knn <- 5

PHATE_data <- phate(data, knn = knn)
PHATE_data <- as_tibble(PHATE_data$embedding)

names(PHATE_data) <- c("emb1", "emb2")

PHATE_data <- PHATE_data |>
  mutate(ID = row_number())

write_rds(PHATE_data, file = paste0("data/limb_muscles/facs_limb_muscles_phate_knn_", knn, ".rds"))


## TriMAP

trimap <- reticulate::import("trimap")

data_vector <- unlist(data)
# Convert the vector into a matrix
data_matrix <- matrix(data_vector, ncol = NCOL(data))

n_inliers <- as.integer(12)
n_outliers <- as.integer(4)
n_random <- as.integer(3)

# Initialize PaCMAP instance
reducer <- trimap$TRIMAP(n_dims = as.integer(2),
                         n_inliers = n_inliers,
                         n_outliers = n_outliers,
                         n_random = n_random)

# Perform dimensionality Reduction
TriMAP_data <- reducer$fit_transform(data_matrix) |>
  as_tibble()

names(TriMAP_data) <- c("emb1", "emb2")

TriMAP_data <- TriMAP_data |>
  mutate(ID = row_number())

write_rds(TriMAP_data, file = paste0("data/limb_muscles/facs_limb_muscles_trimap_n-inliers_", n_inliers, "_n-outliers_", n_outliers, "_n-random_", n_random, ".rds"))

## PaCMAP

pacmap <- reticulate::import("pacmap")

data_vector <- unlist(data)
# Convert the vector into a matrix
data_matrix <- matrix(data_vector, ncol = NCOL(data))

n_neighbors <- as.integer(10)
MN_ratio <- 0.5
FP_ratio <- as.integer(2)
init <- "random"

# Initialize PaCMAP instance
reducer <- pacmap$PaCMAP(n_components = as.integer(2),
                         n_neighbors = n_neighbors,
                         MN_ratio = MN_ratio,
                         FP_ratio = FP_ratio)


# Perform dimensionality Reduction
PacMAP_data <- reducer$fit_transform(data_matrix, init = init) |>
  as_tibble()

names(PacMAP_data) <- c("emb1", "emb2")

PacMAP_data <- PacMAP_data |>
  mutate(ID = row_number())

write_rds(PacMAP_data, file = paste0("data/limb_muscles/facs_limb_muscles_pacmap_n-neighbors_", n_neighbors,"_init_", init, "_MN-ratio_", MN_ratio, "_FP-ratio_", FP_ratio, ".rds"))

## tSNE (best layout)
perplexity <- 15

tSNE_fit <- data |>
  dplyr::select(where(is.numeric)) |>
  Rtsne::Rtsne(perplexity = perplexity,
               pca = FALSE)

tSNE_data <- tSNE_fit$Y |>
  tibble::as_tibble(.name_repair = "unique")
names(tSNE_data) <- c("emb1", "emb2")

tSNE_data <- tSNE_data |>
  mutate(ID = row_number())

write_rds(tSNE_data, file = paste0("data/limb_muscles/facs_limb_muscles_tsne_perplexity_", perplexity, ".rds"))
