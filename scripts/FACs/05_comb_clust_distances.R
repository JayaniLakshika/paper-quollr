library(tidyverse)

dist_df <- read_rds("data/limb_muscles/facs_limb_dist_df.rds")
cluster_df <- read_rds("data/limb_muscles/facs_limb_muscles_cluster_df.rds")

## Map cluster labels for from and to
dist_df <- dplyr::inner_join(dist_df, cluster_df, by = c("from" = "ID")) |>
  dplyr::rename("from_cluster" = "cluster.ids")

dist_df <- dplyr::inner_join(dist_df, cluster_df, by = c("to" = "ID")) |>
  dplyr::rename("to_cluster" = "cluster.ids")

### Run only once
write_rds(dist_df, "data/limb_muscles/facs_limb_dist_df_with clusters.rds")
