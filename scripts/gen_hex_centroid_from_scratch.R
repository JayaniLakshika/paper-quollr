library(quollr)
library(dplyr)
library(ggplot2)

## To compute the diameter of the hexagon
cell_area <- 1
cell_diameter <- sqrt(2 * cell_area / sqrt(3))

# Define hexagonal grid parameters
hex_size <- cell_diameter/2  # Size of the hexagons

buffer_size <- 0.1  # Specify the buffer size here
buffer_size <- hex_size/2  ## should less than hex_size


num_bins_x <- 5 #Exact bynber of bins alng the x-axis

shape_val <- calculate_effective_shape_value(.data = s_curve_noise_umap,
                                             x = "UMAP1", y = "UMAP2")
shape_val

num_bins_y <- 2 * floor(((num_bins_x - 1) * shape_val)/sqrt(3) + 1.5001)
num_bins_y

x_bounds <- seq(min(s_curve_noise_umap$UMAP1) - buffer_size, max(s_curve_noise_umap$UMAP1) + buffer_size, length.out = 5)
y_bounds <- seq(min(s_curve_noise_umap$UMAP2) - buffer_size, max(s_curve_noise_umap$UMAP2) + buffer_size, length.out = 12)


# # Calculate middle points between consecutive values
# middle_points_x <- (x_bounds[-length(x_bounds)] + x_bounds[-1]) / hex_size
# middle_points_y <- (y_bounds[-length(y_bounds)] + y_bounds[-1]) / hex_size


box_points <- expand.grid(x = x_bounds,
                             y = y_bounds)

## For each x-value find the even y value

# Function to generate even y-values for a given x-value
generate_even_y <- function(data) {
  data |>
    dplyr::arrange(y) |>
    dplyr::mutate(is_even = ifelse(seq_along(y) %% 2 == 0, 1, 0))
}

# For each x-value, generate even y-values
box_points <- box_points |>
  dplyr::arrange(x) |>
  dplyr::group_by(x) |>
  dplyr::group_modify(~ generate_even_y(.x))

# split_by_unique_x <- box_points |>
#   dplyr::arrange(x) |>
#   dplyr::group_by(x, .add = TRUE) |>
#   dplyr::group_split()
#
# is_even_vec1 <- split_by_unique_x[[1]] |>
#   dplyr::arrange(y) |>
#   dplyr::mutate(is_even = ifelse(seq_along(y) %% 2 == 0, 1, 0))
#
# is_even_vec2 <- split_by_unique_x[[2]] |>
#   dplyr::arrange(y) |>
#   dplyr::mutate(is_even = ifelse(seq_along(y) %% 2 == 0, 1, 0))
#
# is_even_vec3 <- split_by_unique_x[[3]] |>
#   dplyr::arrange(y) |>
#   dplyr::mutate(is_even = ifelse(seq_along(y) %% 2 == 0, 1, 0))
#
# is_even_vec4 <- split_by_unique_x[[4]] |>
#   dplyr::arrange(y) |>
#   dplyr::mutate(is_even = ifelse(seq_along(y) %% 2 == 0, 1, 0))
#
# is_even_vec5 <- split_by_unique_x[[5]] |>
#   dplyr::arrange(y) |>
#   dplyr::mutate(is_even = ifelse(seq_along(y) %% 2 == 0, 1, 0))
#
#
#
# #box_points$x <- box_points$x + hex_size * ifelse(is_even_vec == 1, 1, 0)
# box_points <- bind_rows(is_even_vec1, is_even_vec2, is_even_vec3, is_even_vec4, is_even_vec5)
x_bounds_n <- box_points$x |> unique()
x_shift <- x_bounds_n[-length(x_bounds_n)] - x_bounds_n[-1]

box_points$x <- box_points$x + x_shift/2 * ifelse(box_points$is_even == 1, 0, 1)

# box_points <- box_points |>
#     dplyr::mutate(is_odd = ifelse(seq_along(box_points$y) %% 2 == 0, 1, 0)) |>
#     dplyr::filter(is_odd == 1)

# middle_box_points <- expand.grid(x = middle_points_x,
#                           y = middle_points_y)
#
# ## Remove all even y rows
#
# box_points <- box_points |>
#   dplyr::mutate(is_odd = ifelse(seq_along(box_points$y) %% 2 == 0, 1, 0)) |>
#   dplyr::filter(is_odd == 0)

# ## Add buffer
# box_points$x_bounds <- box_points$x_bounds + hex_size * sqrt(3)
#
# # Shift every other row of hexagons
# box_points$y_bounds <- box_points$y_bounds + hex_size * ifelse(seq_along(box_points$x) %% 2 == 0, 1, 0)

#box_points$y_bounds <- box_points$y_bounds + hex_size


ggplot() +
  geom_point(data = box_points, aes(x = x, y = y), color = "red")



ggplot() +
  geom_point(data = box_points,
             aes(x = x, y = y, colour = as.factor(is_even)))

# ggplot() +
#   geom_point(data = bind_rows(box_points, middle_box_points), aes(x = x, y = y), color = "red")

#### Full grid


# sx <- 5/diff(x_bounds)
# ## Compute horizontal width of the hexagon
# dx <- 0.5/sx
#
# ## Compute vertical width of the hexagon
# # Adjust for difference in width and height of regular hexagon. 1.15 adjusts
# # for the effect of the overlapping range in y-direction on the resolution
# sy <- (5* shape_val)/diff(y_bounds)
# dy <- (1/sqrt(3))/(2*sy)
#
# ## Obtain hexagon polygon coordinates
# hexC <- hexbin::hexcoords(dx[1], dy[1], n = 1) ## n: Number of hexagons repeat
#
# ## Obtain the number of hexagons in the full grid
# n <- length(box_points$x)
#
# ## Generate the size vector of the hexagons (since regular hexagons)
# size <- rep(1, length(box_points$x))
#
# ## Generate the coordinates for the hexagons
# full_hex_coords <- tibble::tibble( x = rep.int(hexC$x, n) * rep(size, each = 6) + rep(box_points$x, each = 6),
#                                    y = rep.int(hexC$y, n) * rep(size, each = 6) + rep(box_points$y, each = 6),
#                                    id = rep(1:length(box_points$x), each = 6))
#







#box_points <- bind_rows(box_points, middle_box_points)
## Generate all coordinates of hexagons
hex_grid_n <- full_hex_grid(box_points)

ggplot(data = hex_grid_n, aes(x = x, y = y)) + geom_polygon(fill = "white", color = "black", aes(group = id)) +
  geom_point(data = box_points, aes(x = x, y = y), color = "red") +
  geom_point(data = s_curve_noise_umap, aes(x = UMAP1, y = UMAP2), color = "blue")

### Assign data to centroids

find_nearest_centroid <- function(data, centroids) {
  distances <- as.matrix(stats::dist(rbind(data$UMAP1, data$UMAP2), rbind(centroids$x, centroids$y), method = "euclidean"))
  nearest_centroids <- apply(distances, 1, which.min)
  return(nearest_centroids)
}

find_nearest_centroid(s_curve_noise_umap, box_points)

### Map hexIDs

full_centroid_df <- hex_grid_n

vec1 <- stats::setNames(rep("", 2), c("x", "y"))  ## Define column names

## Define a dataset to store all the centroids with the respective coordinates
full_grid_with_hexbin_id <- dplyr::bind_rows(vec1)[0, ]
full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
  dplyr::mutate_if(is.character, as.numeric)

for(i in 1:length(sort(unique(full_centroid_df$y)))){

  ## Filter the data set with specific y value
  specific_y_val_df <- full_centroid_df |>
    dplyr::filter(y == sort(unique(full_centroid_df$y))[i])

  ## orderd the x values
  ordered_x_df <- specific_y_val_df |>
    dplyr::arrange(x)

  full_grid_with_hexbin_id <- dplyr::bind_rows(full_grid_with_hexbin_id, ordered_x_df)

}

## Add the column with hexagonal bin ID
full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
  dplyr::mutate(hexID = dplyr::row_number())

## Rename columns
full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
  dplyr::rename("c_x" = "x",
                "c_y" = "y")

ggplot() +
  geom_text(data = full_grid_with_hexbin_id, aes(x = c_x, y = c_y, label = hexID), color = "red")


hex_full_count_df <- inner_join(full_grid_with_hexbin_id, hex_grid_n)

ggplot(data = hex_full_count_df, aes(x = x, y = y)) +
  geom_polygon(color = "black", aes(group = id), colour = "white") +
  geom_point(data = box_points, aes(x = x, y = y), colour = "red")

### Assign points
first_polygon <- hex_grid_n |> dplyr::filter(id == 21)
range(first_polygon$x)
range(first_polygon$y)

filtered_data_points <- s_curve_noise_umap |>
  dplyr::filter((UMAP1 >= 0.5) & (UMAP1 <= 1.8)) |>
  dplyr::filter((UMAP2 >= 0.3) & (UMAP2 <= 1))

### First find minimum and max values of UMAP data
### Find the box of centroids within that values (start the search from there)
### then check of one by one hexagons by taking their min and max polygon coordinates
### Since min and max values can be take some part of nearest y hexagons, need to search on up and down 4 neighbours
### by computing distances to the centroids
### Filter the UMAP data for selected min and max polygon coordinates, and start the search for that hexgaon and their neighbors
### Assign the hexID with minimum distance
