library(readr)
library(quollr)
library(dplyr)

set.seed(20240110)

data <- read_rds("data/limb_muscles/facs_limb_muscles_pcs_10.rds")

## For umap
umap_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_umap_n-neigbors_15_min-dist_0.1.rds")

error_limb_muscles_umap <- gen_diffbin1_errors(highd_data = data, nldr_data = umap_limb_muscles) |>
  dplyr::mutate(method = "UMAP_15_min_dist_0.1")

write_rds(error_limb_muscles_umap, "data/limb_muscles/error_limb_muscles_umap_n-neigbors_15_min-dist_0.1.rds")

###########

## For tsne
tsne_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_tsne_perplexity_30.rds")

error_limb_muscles_tsne <- gen_diffbin1_errors(highd_data = data, nldr_data = tsne_limb_muscles) |>
  dplyr::mutate(method = "tsne_30")

write_rds(error_limb_muscles_tsne, "data/limb_muscles/error_limb_muscles_tsne_perplexity_30.rds")

###########

## For phate
phate_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_phate_knn_5.rds")

error_limb_muscles_phate <- gen_diffbin1_errors(highd_data = data, nldr_data = phate_limb_muscles) |>
  dplyr::mutate(method = "phate_5")

write_rds(error_limb_muscles_phate, "data/limb_muscles/error_limb_muscles_phate_knn_5.rds")

###########

## For trimap
trimap_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_trimap_n-inliers_12_n-outliers_4_n-random_3.rds")

error_limb_muscles_trimap <- gen_diffbin1_errors(highd_data = data, nldr_data = trimap_limb_muscles) |>
  dplyr::mutate(method = "trimap_n-inliers_12_n-outliers_4_n-random_3")

write_rds(error_limb_muscles_trimap, "data/limb_muscles/error_limb_muscles_trimap_n-inliers_12_n-outliers_4_n-random_3.rds")

###########

## For pacmap
pacmap_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_pacmap_n-neighbors_10_init_random_MN-ratio_0.5_FP-ratio_2.rds")

error_limb_muscles_pacmap <- gen_diffbin1_errors(highd_data = data, nldr_data = pacmap_limb_muscles) |>
  dplyr::mutate(method = "pacmap_n-neighbors_10_init_random_MN-ratio_0.5_FP-ratio_2")

write_rds(error_limb_muscles_pacmap, "data/limb_muscles/error_limb_muscles_pacmap_n-neighbors_10_init_random_MN-ratio_0.5_FP-ratio_2.rds")

###########

## For trimap
trimap_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_trimap_n-inliers_5_n-outliers_4_n-random_3.rds")

error_limb_muscles_trimap <- gen_diffbin1_errors(highd_data = data, nldr_data = trimap_limb_muscles) |>
  dplyr::mutate(method = "trimap_n-inliers_5_n-outliers_4_n-random_3")

write_rds(error_limb_muscles_trimap, "data/limb_muscles/error_limb_muscles_trimap_n-inliers_5_n-outliers_4_n-random_3.rds")

## For umap
umap_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_umap_n-neigbors_30_min-dist_0.1.rds")

error_limb_muscles_umap <- gen_diffbin1_errors(highd_data = data, nldr_data = umap_limb_muscles) |>
  dplyr::mutate(method = "UMAP_30_min_dist_0.13")

write_rds(error_limb_muscles_umap, "data/limb_muscles/error_limb_muscles_umap_n-neigbors_30_min-dist_0.1.rds")

###########

## For tsne
tsne_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_tsne_perplexity_15.rds")

error_limb_muscles_tsne <- gen_diffbin1_errors(highd_data = data, nldr_data = tsne_limb_muscles) |>
  dplyr::mutate(method = "tsne_15")

write_rds(error_limb_muscles_tsne, "data/limb_muscles/error_limb_muscles_tsne_perplexity_15.rds")
