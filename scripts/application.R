<!--
  ### Appropriate hex size value to have regular hexagons

  ```{r}

min_value_x <- min(s_curve_noise_umap$UMAP1)
max_value_x <- max(s_curve_noise_umap$UMAP1)
min_value_y <- min(s_curve_noise_umap$UMAP2)
max_value_y <- max(s_curve_noise_umap$UMAP2)

# Calculate the maximum span along x and y axes
x_span <- max_value_x - min_value_x
y_span <- max_value_y - min_value_y

hex_size <- min(x_span, y_span)/5.5

num_bins_x <- calculate_effective_x_bins(.data = s_curve_noise_umap, x = "UMAP1", hex_size = hex_size)
num_bins_y <- calculate_effective_y_bins(.data = s_curve_noise_umap, y = "UMAP2", hex_size = hex_size)

all_centroids_df <- generate_full_grid_centroids(nldr_df = s_curve_noise_umap,
                                                 x = "UMAP1", y = "UMAP2",
                                                 num_bins_x = num_bins_x,
                                                 num_bins_y = num_bins_y,
                                                 buffer_size = NA, hex_size = hex_size)
hex_grid <- gen_hex_coordinates(all_centroids_df, hex_size = hex_size)
ggplot(data = hex_grid, aes(x = x, y = y)) + geom_polygon(fill = "white", color = "black", aes(group = id)) +
  geom_point(data = all_centroids_df, aes(x = x, y = y), color = "red")
```
-->


### Model function

`fit_high_d_model()` function is used to generate the 2D model and the high-D model.

```{r, eval=FALSE}
fit_high_d_model(training_data = training_data, nldr_df_with_id = s_curve_noise_umap, x = "UMAP1",
                 y = "UMAP1", num_bins_x = NA, num_bins_y = NA,
                 hex_size = NA, buffer_size = NA,
                 is_bin_centroid = TRUE,
                 is_rm_lwd_hex = FALSE,
                 benchmark_to_rm_lwd_hex = NA,
                 is_avg_high_d = TRUE, column_start_text = "x")
```

#### Benchmark value to remove the low-density hexagons

```{r, eval=FALSE}
## As an option first quantile considered as a default
benchmark_to_rm_lwd_hex <- quantile(df_bin_centroids$std_counts)[2] + 0.01

## To identify low density hexagons
df_bin_centroids_low <- df_bin_centroids |>
  dplyr::filter(std_counts <= benchmark_to_rm_lwd_hex)

## To identify low-density hexagons needed to remove by investigating neighbouring mean density
identify_rm_bins <- find_low_density_hexagons(df_bin_centroids_all = df_bin_centroids, num_bins_x = num_bins_x,
                                              df_bin_centroids_low = df_bin_centroids_low)
```


## Shift

```{r}
all_centroids_df_shift <- extract_coord_of_shifted_hex_grid(nldr_df = s_curve_noise_umap,
                                                            x = "UMAP1", y = "UMAP2",
                                                            num_bins_x = num_bins_x,
                                                            num_bins_y = num_bins_y,
                                                            shift_x = 0.2690002, shift_y = 0.271183,
                                                            buffer_size = NA, hex_size = NA)

glimpse(all_centroids_df_shift)
```

```{r}
hex_grid <- gen_hex_coordinates(all_centroids_df_shift)
glimpse(hex_grid)
```

```{r}
ggplot(data = hex_grid, aes(x = x, y = y)) + geom_polygon(fill = "white", color = "black", aes(group = id)) +
  geom_point(data = all_centroids_df_shift, aes(x = x, y = y), color = "red")
```

```{r}
full_grid_with_hexbin_id <- map_hexbin_id(all_centroids_df_shift)

ggplot(data = hex_grid, aes(x = x, y = y)) + geom_polygon(fill = "white", color = "black", aes(group = id)) +
  geom_text(data = full_grid_with_hexbin_id, aes(x = c_x, y = c_y, label = hexID))
```

```{r}
full_grid_with_polygon_id <- map_polygon_id(full_grid_with_hexbin_id, hex_grid)
```

```{r}
s_curve_noise_umap_with_id <- assign_data(s_curve_noise_umap, full_grid_with_hexbin_id)
```

```{r}
df_with_std_counts <- compute_std_counts(nldr_df = s_curve_noise_umap_with_id)
```

```{r}
hex_full_count_df <- generate_full_grid_info(full_grid_with_polygon_id, df_with_std_counts, hex_grid)
```

```{r}
ggplot(data = hex_grid, aes(x = x, y = y)) + geom_polygon(fill = "white", color = "black", aes(group = id)) +
  geom_point(data = s_curve_noise_umap, aes(x = UMAP1, y = UMAP2), color = "blue")
```

```{r}
ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_text(aes(x = c_x, y = c_y, label = hexID)) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff")

```

```{r}
df_bin_centroids <- hex_full_count_df[complete.cases(hex_full_count_df[["std_counts"]]), ] |>
  dplyr::select("c_x", "c_y", "hexID", "std_counts") |>
  dplyr::distinct() |>
  dplyr::rename(c("x" = "c_x", "y" = "c_y"))

df_bin_centroids
```

```{r}
tr1_object <- triangulate_bin_centroids(df_bin_centroids, x = "x", y = "y")
tr_from_to_df <- generate_edge_info(triangular_object = tr1_object)
```

```{r, eval=FALSE}
bin_centroids_shift <- ggplot(data = hex_full_count_df, aes(x = c_x, y = c_y)) +
  geom_point(color = "#bdbdbd") +
  geom_point(data = shifted_hex_coord_df, aes(x = c_x, y = c_y), color = "#feb24c") +
  coord_cartesian(xlim = c(-5, 8), ylim = c(-10, 10)) +
  theme_void() +
  theme(legend.position="none", legend.direction="horizontal", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6)) +
  guides(fill = guide_colourbar(title = "Standardized count")) +
  annotate(geom = 'text', label = "a", x = -Inf, y = Inf, hjust = -0.3, vjust = 1, size = 3)

hex_grid_shift <- ggplot(data = shifted_hex_coord_df, aes(x = x, y = y)) +
  geom_polygon(fill = NA, color = "#feb24c", aes(group = polygon_id)) +
  geom_polygon(data = hex_full_count_df, aes(x = x, y = y, group = polygon_id),
               fill = NA, color = "#bdbdbd") +
  coord_cartesian(xlim = c(-5, 8), ylim = c(-10, 10)) +
  theme_void() +
  theme(legend.position="none", legend.direction="horizontal", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6)) +
  guides(fill = guide_colourbar(title = "Standardized count")) +
  annotate(geom = 'text', label = "b", x = -Inf, y = Inf, hjust = -0.3, vjust = 1, size = 3)

## Before shift
before_shift_plot <- ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_text(aes(x = c_x, y = c_y, label = hexID), size = 2) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff", option = "C") +
  coord_equal() +
  theme_void() +
  theme(legend.position="bottom", legend.direction="horizontal", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6)) +
  guides(fill = guide_colourbar(title = "Standardized count")) +
  annotate(geom = 'text', label = "a", x = -Inf, y = Inf, hjust = -0.3, vjust = 1, size = 3)


## After shift
after_shift_plot <- ggplot(data = shifted_hex_coord_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_text(aes(x = c_x, y = c_y, label = hexID), size = 2) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff", option = "C") +
  coord_equal() +
  theme_void() +
  theme(legend.position="none", legend.direction="horizontal", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6)) +
  guides(fill = guide_colourbar(title = "Standardized count")) +
  annotate(geom = 'text', label = "b", x = -Inf, y = Inf, hjust = -0.3, vjust = 1, size = 3)

```

## Application
```{r}
medlea_df <- read_csv("data/medlea_dataset.csv")
names(medlea_df)[2:(NCOL(medlea_df) - 1)] <- paste0("x", 1:(NCOL(medlea_df) - 2))

medlea_df <- medlea_df |> ## Since only contains zeros
  select(-x10)

#medlea_df[,2:(NCOL(medlea_df) - 1)] <- scale(medlea_df[,2:(NCOL(medlea_df) - 1)])

calculate_pca <- function(feature_dataset, num_pcs){
  pcaY_cal <- prcomp(feature_dataset, center = TRUE, scale = TRUE)
  PCAresults <- data.frame(pcaY_cal$x[, 1:num_pcs])
  summary_pca <- summary(pcaY_cal)
  var_explained_df <- data.frame(PC= paste0("PC",1:50),
                                 var_explained=(pcaY_cal$sdev[1:50])^2/sum((pcaY_cal$sdev[1:50])^2))
  return(list(prcomp_out = pcaY_cal,pca_components = PCAresults, summary = summary_pca, var_explained_pca  = var_explained_df))
}
features <- medlea_df[,2:(NCOL(medlea_df) - 1)]
pca_ref_calc <- calculate_pca(features, 8)
pca_ref_calc$summary

var_explained_df <- pca_ref_calc$var_explained_pca
data_pca <- pca_ref_calc$pca_components |>
  mutate(ID = 1:NROW(pca_ref_calc$pca_components),
         shape_label = medlea_df$Shape_label)

var_explained_df |>
  ggplot(aes(x = PC,y = var_explained, group = 1))+
  geom_point(size=1)+
  geom_line()+
  labs(title="Scree plot: PCA on scaled data") +
  scale_x_discrete(limits = paste0(rep("PC", 50), 1:50)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

data_split <- initial_split(data_pca)
training_data <- training(data_split) |>
  arrange(ID)
test_data <- testing(data_split) |>
  arrange(ID)
```

```{r}
UMAP_fit <- umap(training_data |> dplyr::select(-c(ID, shape_label)), n_neighbors = 37, n_components =  2)

UMAP_data <- UMAP_fit$layout |>
  as.data.frame()
names(UMAP_data)[1:(ncol(UMAP_data))] <- paste0(rep("UMAP",(ncol(UMAP_data))), 1:(ncol(UMAP_data)))

UMAP_data <- UMAP_data |>
  mutate(ID = training_data$id)

UMAP_data_with_label <- UMAP_data |>
  mutate(shape_label = training_data$shape_label)
```

```{r}
UMAP_data_with_label |>
  ggplot(aes(x = UMAP1,
             y = UMAP2, color = shape_label))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + #ggtitle("(a)") +
  theme_linedraw() +
  theme(legend.position = "none", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm')) +
  scale_color_manual(values=c("#b15928", "#1f78b4", "#cab2d6", "#ccebc5", "#fb9a99", "#e31a1c", "#6a3d9a", "#ff7f00", "#ffed6f", "#fdbf6f", "#ffff99", "#a6cee3", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#b2df8a", "#bc80bd", "#33a02c", "#ccebc5", "#ffed6f", "#000000", "#bdbdbd"))
```

```{r}
tSNE_data <- Fit_tSNE(training_data |> dplyr::select(-c(ID, shape_label)), opt_perplexity = calculate_effective_perplexity(training_data |> dplyr::select(-c(ID, shape_label))), with_seed = 20240110)

tSNE_data <- tSNE_data |>
  select(-ID) |>
  mutate(ID = training_data$ID)

tSNE_data_with_label <- tSNE_data |>
  mutate(shape_label = training_data$shape_label)

tSNE_data_with_label |>
  ggplot(aes(x = tSNE1,
             y = tSNE2, color = shape_label))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + #ggtitle("(a)") +
  theme_linedraw() +
  theme(legend.position = "none", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm')) +
  scale_color_manual(values=c("#b15928", "#1f78b4", "#cab2d6", "#ccebc5", "#fb9a99", "#e31a1c", "#6a3d9a", "#ff7f00", "#ffed6f", "#fdbf6f", "#ffff99", "#a6cee3", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#b2df8a", "#bc80bd", "#33a02c", "#ccebc5", "#ffed6f", "#000000", "#bdbdbd"))
```

```{r}
PHATE_data <- Fit_PHATE(training_data |> dplyr::select(-c(ID, shape_label)), knn = 5, with_seed = 20240110)
PHATE_data <- PHATE_data |>
  select(PHATE1, PHATE2)
PHATE_data <- PHATE_data |>
  mutate(ID = training_data$ID)

PHATE_data_with_label <- PHATE_data |>
  mutate(shape_label = training_data$shape_label)

PHATE_data_with_label |>
  ggplot(aes(x = PHATE1,
             y = PHATE2, color = shape_label))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + #ggtitle("(a)") +
  theme_linedraw() +
  theme(legend.position = "none", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm')) +
  scale_color_manual(values=c("#b15928", "#1f78b4", "#cab2d6", "#ccebc5", "#fb9a99", "#e31a1c", "#6a3d9a", "#ff7f00", "#ffed6f", "#fdbf6f", "#ffff99", "#a6cee3", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#b2df8a", "#bc80bd", "#33a02c", "#ccebc5", "#ffed6f", "#000000", "#bdbdbd"))
```

```{r}
tem_dir <- tempdir()

Fit_TriMAP_data(training_data |> dplyr::select(-c(ID, shape_label)), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_TriMAP_values.csv")

Fit_TriMAP(as.integer(2), as.integer(5), as.integer(4), as.integer(3), path, path2)

TriMAP_data <- read_csv(path2)
TriMAP_data <- TriMAP_data |>
  mutate(ID = training_data$ID)

TriMAP_data_with_label <- TriMAP_data |>
  mutate(shape_label = training_data$shape_label)

TriMAP_data_with_label |>
  ggplot(aes(x = TriMAP1,
             y = TriMAP2, color = shape_label))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + #ggtitle("(a)") +
  theme_linedraw() +
  theme(legend.position = "none", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm')) +
  scale_color_manual(values=c("#b15928", "#1f78b4", "#cab2d6", "#ccebc5", "#fb9a99", "#e31a1c", "#6a3d9a", "#ff7f00", "#ffed6f", "#fdbf6f", "#ffff99", "#a6cee3", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#b2df8a", "#bc80bd", "#33a02c", "#ccebc5", "#ffed6f", "#000000", "#bdbdbd"))

```


```{r}
tem_dir <- tempdir()

Fit_PacMAP_data(training_data |> dplyr::select(-c(ID, shape_label)), tem_dir)

path <- file.path(tem_dir, "df_2_without_class.csv")
path2 <- file.path(tem_dir, "dataset_3_PaCMAP_values.csv")

Fit_PaCMAP(as.integer(2), as.integer(10), "random", 0.9, as.integer(2), path, path2)

PaCMAP_data <- read_csv(path2)
PaCMAP_data <- PaCMAP_data |>
  mutate(ID = training_data$ID)

PaCMAP_data_with_label <- PaCMAP_data |>
  mutate(shape_label = training_data$shape_label)

PaCMAP_data_with_label |>
  ggplot(aes(x = PaCMAP1,
             y = PaCMAP2, color = shape_label))+
  geom_point(alpha=0.5) +
  coord_equal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + #ggtitle("(a)") +
  theme_linedraw() +
  theme(legend.position = "none", plot.title = element_text(size = 7, hjust = 0.5, vjust = -0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #change legend key width
        legend.title = element_text(size=5), #change legend title font size
        legend.text = element_text(size=4),
        legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm')) +
  scale_color_manual(values=c("#b15928", "#1f78b4", "#cab2d6", "#ccebc5", "#fb9a99", "#e31a1c", "#6a3d9a", "#ff7f00", "#ffed6f", "#fdbf6f", "#ffff99", "#a6cee3", "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#b2df8a", "#bc80bd", "#33a02c", "#ccebc5", "#ffed6f", "#000000", "#bdbdbd"))
```

```{r}
num_bins_x <- calculate_effective_x_bins(.data = tSNE_data, x = "tSNE1", hex_size = NA)
num_bins_x
```

```{r}
num_bins_y <- calculate_effective_y_bins(.data = tSNE_data, y = "tSNE2", hex_size = NA)
num_bins_y
```

```{r}
all_centroids_df <- generate_full_grid_centroids(nldr_df = tSNE_data,
                                                 x = "tSNE1", y = "tSNE2",
                                                 num_bins_x = num_bins_x,
                                                 num_bins_y = num_bins_y,
                                                 buffer_size = NA, hex_size = NA)


hex_grid <- gen_hex_coordinates(all_centroids_df)

full_grid_with_hexbin_id <- map_hexbin_id(all_centroids_df)

full_grid_with_polygon_id <- map_polygon_id(full_grid_with_hexbin_id, hex_grid)

tSNE_data_with_id <- assign_data(tSNE_data, full_grid_with_hexbin_id)

df_with_std_counts <- compute_std_counts(nldr_df = tSNE_data_with_id)

hex_full_count_df <- generate_full_grid_info(full_grid_with_polygon_id, df_with_std_counts, hex_grid)

ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = polygon_id, fill = std_counts)) +
  geom_text(aes(x = c_x, y = c_y, label = hexID), size = 2) +
  scale_fill_viridis_c(direction = -1, na.value = "#ffffff")
```


```{r}
ggplot(data = hex_grid, aes(x = x, y = y)) + geom_polygon(fill = "white", color = "black", aes(group = id)) +
  geom_point(data = tSNE_data, aes(x = tSNE1, y = tSNE2), color = "blue")
```

```{r}
df_bin_centroids <- extract_hexbin_centroids(hex_full_count_df)
```

```{r}
tr1_object <- triangulate_bin_centroids(df_bin_centroids, x = "x", y = "y")
tr_from_to_df <- generate_edge_info(triangular_object = tr1_object)
```

```{r}
## To generate a data set with high-D and 2D training data
df_all <- training_data |> dplyr::select(-c(ID, shape_label)) |>
  dplyr::bind_cols(tSNE_data_with_id)

## To generate averaged high-D data

df_bin <- avg_highD_data(.data = df_all, column_start_text = "PC") ## Need to pass ID column name
```

```{r}
## Compute 2D distances
distance <- cal_2d_dist(tr_from_to_df_coord = tr_from_to_df)

plot_dist(distance)

benchmark <- find_benchmark_value(distance_edges = distance, distance_col = "distance")
```

```{r}
trimesh <- ggplot(df_bin_centroids, aes(x = x, y = y)) +
  geom_point(size = 0.1) +
  geom_trimesh() +
  coord_equal()

trimesh
```

```{r}
trimesh_gr <- colour_long_edges(distance_edges = distance, benchmark_value = benchmark,
                                tr_from_to_df_coord = tr_from_to_df, distance_col = "distance")

trimesh_gr
```

```{r}
trimesh_removed <- remove_long_edges(distance_edges = distance, benchmark_value = benchmark,
                                     tr_from_to_df_coord = tr_from_to_df, distance_col = "distance")
trimesh_removed
```

```{r}
tour1 <- show_langevitour(df_all, df_bin, df_bin_centroids, benchmark_value = benchmark,
                          distance = distance, distance_col = "distance", column_start_text = "PC")
tour1
```
