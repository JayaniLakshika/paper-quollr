library(readr)
library(quollr)
library(dplyr)

data <- read_rds("data/limb_muscles/facs_limb_muscles_pcs_10.rds")

## For umap
umap_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_umap_n-neigbors_15_min-dist_0.1.rds")

limb_scaled_obj_umap <- gen_scaled_data(
  data = umap_limb_muscles)
umap_limb_muscles_scaled <- limb_scaled_obj_umap$scaled_nldr

lim1 <- limb_scaled_obj_umap$lim1
lim2 <- limb_scaled_obj_umap$lim2
r2_umap <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_limb_muscles <- 2:37 #sqrt(NROW(data)/r2_umap)

error_limb_muscles_umap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_limb_muscles) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_umap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  limb_model <- fit_highd_model(
    highd_data = data,
    nldr_data = umap_limb_muscles_scaled,
    bin1 = xbins,
    r2 = r2_umap,
    q = 0.1,
    is_bin_centroid = TRUE
  )

  df_bin_centroids_limb_muscles <- limb_model$df_bin_centroids
  df_bin_limb_muscles <- limb_model$df_bin

  ## Compute error
  error_df <- glance(
    model_2d = df_bin_centroids_limb_muscles,
    model_highd = df_bin_limb_muscles,
    highd_data = data) |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_limb_muscles),
           method = "UMAP_30_min_dist_0.1",
           a1 = a1)

  error_limb_muscles_umap <- bind_rows(error_limb_muscles_umap, error_df)

}

write_rds(error_limb_muscles_umap, "data/limb_muscles/error_limb_muscles_umap_n-neigbors_30_min-dist_0.1.rds")

###########

## For tsne
tsne_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_tsne_perplexity_30.rds")

limb_scaled_obj_tsne <- gen_scaled_data(
  data = tsne_limb_muscles)
tsne_limb_muscles_scaled <- limb_scaled_obj_tsne$scaled_nldr

lim1 <- limb_scaled_obj_tsne$lim1
lim2 <- limb_scaled_obj_tsne$lim2
r2_tsne <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_limb_muscles <- 2:34 #sqrt(NROW(data)/r2_tsne)

error_limb_muscles_tsne <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_limb_muscles) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_tsne, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  limb_model <- fit_highd_model(
    highd_data = data,
    nldr_data = tsne_limb_muscles_scaled,
    bin1 = xbins,
    r2 = r2_tsne,
    q = 0.1,
    is_bin_centroid = TRUE
  )

  df_bin_centroids_limb_muscles <- limb_model$df_bin_centroids
  df_bin_limb_muscles <- limb_model$df_bin

  ## Compute error
  error_df <- glance(
    model_2d = df_bin_centroids_limb_muscles,
    model_highd = df_bin_limb_muscles,
    highd_data = data) |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_limb_muscles),
           method = "tsne_30",
           a1 = a1)

  error_limb_muscles_tsne <- bind_rows(error_limb_muscles_tsne, error_df)

}

write_rds(error_limb_muscles_tsne, "data/limb_muscles/error_limb_muscles_tsne_perplexity_30.rds")

###########

## For phate
phate_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_phate_knn_5.rds")

limb_scaled_obj_phate <- gen_scaled_data(
  data = phate_limb_muscles)
phate_limb_muscles_scaled <- limb_scaled_obj_phate$scaled_nldr

lim1 <- limb_scaled_obj_phate$lim1
lim2 <- limb_scaled_obj_phate$lim2
r2_phate <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_limb_muscles <- 2:34 #sqrt(NROW(data)/r2_phate)

error_limb_muscles_phate <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_limb_muscles) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_phate, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  limb_model <- fit_highd_model(
    highd_data = data,
    nldr_data = phate_limb_muscles_scaled,
    bin1 = xbins,
    r2 = r2_phate,
    q = 0.1,
    is_bin_centroid = TRUE
  )

  df_bin_centroids_limb_muscles <- limb_model$df_bin_centroids
  df_bin_limb_muscles <- limb_model$df_bin

  ## Compute error
  error_df <- glance(
    model_2d = df_bin_centroids_limb_muscles,
    model_highd = df_bin_limb_muscles,
    highd_data = data) |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_limb_muscles),
           method = "phate_5",
           a1 = a1)

  error_limb_muscles_phate <- bind_rows(error_limb_muscles_phate, error_df)

}

write_rds(error_limb_muscles_phate, "data/limb_muscles/error_limb_muscles_phate_knn_5.rds")

###########

## For trimap
trimap_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_trimap_n-inliers_12_n-outliers_4_n-random_3.rds")

limb_scaled_obj_trimap <- gen_scaled_data(
  data = trimap_limb_muscles)
trimap_limb_muscles_scaled <- limb_scaled_obj_trimap$scaled_nldr

lim1 <- limb_scaled_obj_trimap$lim1
lim2 <- limb_scaled_obj_trimap$lim2
r2_trimap <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_limb_muscles <- 2:129 #sqrt(NROW(data)/r2_trimap)

error_limb_muscles_trimap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_limb_muscles) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_trimap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  limb_model <- fit_highd_model(
    highd_data = data,
    nldr_data = trimap_limb_muscles_scaled,
    bin1 = xbins,
    r2 = r2_trimap,
    q = 0.1,
    is_bin_centroid = TRUE
  )

  df_bin_centroids_limb_muscles <- limb_model$df_bin_centroids
  df_bin_limb_muscles <- limb_model$df_bin

  ## Compute error
  error_df <- glance(
    model_2d = df_bin_centroids_limb_muscles,
    model_highd = df_bin_limb_muscles,
    highd_data = data) |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_limb_muscles),
           method = "trimap_n-inliers_12_n-outliers_4_n-random_3",
           a1 = a1)

  error_limb_muscles_trimap <- bind_rows(error_limb_muscles_trimap, error_df)

}

write_rds(error_limb_muscles_trimap, "data/limb_muscles/error_limb_muscles_trimap_n-inliers_12_n-outliers_4_n-random_3.rds")

###########

## For pacmap
pacmap_limb_muscles <- read_rds("data/limb_muscles/facs_limb_muscles_pacmap_n-neighbors_10_init_random_MN-ratio_0.5_FP-ratio_2.rds")

limb_scaled_obj_pacmap <- gen_scaled_data(
  data = pacmap_limb_muscles)
pacmap_limb_muscles_scaled <- limb_scaled_obj_pacmap$scaled_nldr

lim1 <- limb_scaled_obj_pacmap$lim1
lim2 <- limb_scaled_obj_pacmap$lim2
r2_pacmap <- diff(lim2)/diff(lim1)

## To initialize number of bins along the x-axis
bin1_vec_limb_muscles <- 2:38 #sqrt(NROW(data)/r2_pacmap)

error_limb_muscles_pacmap <- data.frame(matrix(nrow = 0, ncol = 0))

for (xbins in bin1_vec_limb_muscles) {

  hb_obj <- calc_bins_y(bin1 = xbins, r2 = r2_pacmap, q = 0.1)

  bin2 <- hb_obj$bin2
  a1 <- hb_obj$a1

  limb_model <- fit_highd_model(
    highd_data = data,
    nldr_data = pacmap_limb_muscles_scaled,
    bin1 = xbins,
    r2 = r2_pacmap,
    q = 0.1,
    is_bin_centroid = TRUE
  )

  df_bin_centroids_limb_muscles <- limb_model$df_bin_centroids
  df_bin_limb_muscles <- limb_model$df_bin

  ## Compute error
  error_df <- glance(
    model_2d = df_bin_centroids_limb_muscles,
    model_highd = df_bin_limb_muscles,
    highd_data = data) |>
    mutate(bin1 = xbins,
           bin2 = bin2,
           b = bin1 * bin2,
           b_non_empty = NROW(df_bin_centroids_limb_muscles),
           method = "pacmap_n-neighbors_10_init_random_MN-ratio_0.5_FP-ratio_2",
           a1 = a1)

  error_limb_muscles_pacmap <- bind_rows(error_limb_muscles_pacmap, error_df)

}

write_rds(error_limb_muscles_pacmap, "data/limb_muscles/error_limb_muscles_pacmap_n-neighbors_10_init_random_MN-ratio_0.5_FP-ratio_2.rds")


