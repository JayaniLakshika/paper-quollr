library(tidyverse)

dist_df <- read_rds("data/limb_muscles/facs_limb_dist_df_with clusters.rds")

## Within clusters

### For cluster 0

dist_df_clust0 <- dist_df |>
  filter(from_cluster == 0) |>
  filter(to_cluster == 0)

dist_df_clust0 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_30)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

dist_df_clust0 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_15)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

### For cluster 1

dist_df_clust1 <- dist_df |>
  filter(from_cluster == 1) |>
  filter(to_cluster == 1)

dist_df_clust1 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_30)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

dist_df_clust1 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_15)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

### For cluster 2

dist_df_clust2 <- dist_df |>
  filter(from_cluster == 2) |>
  filter(to_cluster == 2)

dist_df_clust2 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_30)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

dist_df_clust2 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_15)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

### For cluster 3

dist_df_clust3 <- dist_df |>
  filter(from_cluster == 3) |>
  filter(to_cluster == 3)

dist_df_clust3 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_30)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

dist_df_clust3 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_15)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

### For cluster 4

dist_df_clust4 <- dist_df |>
  filter(from_cluster == 4) |>
  filter(to_cluster == 4)

dist_df_clust4 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_30)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

dist_df_clust4 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_15)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

### For cluster 5

dist_df_clust5 <- dist_df |>
  filter(from_cluster == 5) |>
  filter(to_cluster == 5)

dist_df_clust5 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_30)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

dist_df_clust5 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_15)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

### For cluster 6

dist_df_clust6 <- dist_df |>
  filter(from_cluster == 6) |>
  filter(to_cluster == 6)

dist_df_clust6 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_30)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

dist_df_clust6 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_15)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

## Between clusters

### Between 0 and 2

dist_df_clust02 <- dist_df |>
  filter(from_cluster == 0) |>
  filter(to_cluster == 2)

dist_df_clust02 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_30)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

dist_df_clust02 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_15)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

###
dist_df_clust02 <- dist_df |>
  filter(from_cluster == 2) |>
  filter(to_cluster == 0)

dist_df_clust02 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_30)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

dist_df_clust02 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_15)) +
  geom_point(alpha=0.5) +
  theme(aspect.ratio = 1)

dist_df_clust02 |>
  ggplot(aes(x = dist_highd,
             y = dist_2d_tsne_30)) +
  geom_point(alpha=0.1) +
  geom_point(aes(x = dist_highd,
                 y = dist_2d_tsne_15),
             alpha=0.1,
             color = "#8da0cb") +
  theme(aspect.ratio = 1)

## Add interactivity with 2D layouts, distances comp, highd data
