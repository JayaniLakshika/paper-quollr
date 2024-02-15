library(quollr)
library(dplyr)
library(ggplot2)

## To compute the diameter of the hexagon
cell_area <- 1
cell_diameter <- sqrt(2 * cell_area / sqrt(3))

# Define hexagonal grid parameters
hex_size <- cell_diameter/2  # Size of the hexagons

buffer_size <- 0.01  # Specify the buffer size here
buffer_size <- hex_size/2  ## should less than hex_size


num_bins_x <- 5 #Exact bynber of bins alng the x-axis

shape_val <- calculate_effective_shape_value(.data = s_curve_noise_umap,
                                             x = "UMAP1", y = "UMAP2")
shape_val

num_bins_y <- 2 * floor(((num_bins_x - 1) * shape_val)/sqrt(3) + 1.5001)
num_bins_y

## Maximum shift is cell_diameter
shift_x <- 1.051547
shift_y <- 0.2686425

x_bounds <- seq(min(s_curve_noise_umap$UMAP1) - buffer_size + shift_x, max(s_curve_noise_umap$UMAP1) + buffer_size + shift_x, length.out = 5)
y_bounds <- seq(min(s_curve_noise_umap$UMAP2) - buffer_size + shift_y, max(s_curve_noise_umap$UMAP2) + buffer_size + shift_y, length.out = 12)


# # Calculate middle points between consecutive values
# middle_points_x <- (x_bounds[-length(x_bounds)] + x_bounds[-1]) / hex_size
# middle_points_y <- (y_bounds[-length(y_bounds)] + y_bounds[-1]) / hex_size


box_points <- expand.grid(x = x_bounds,
                          y = y_bounds)

## For each x-value find the even y value

split_by_unique_x <- box_points |>
  dplyr::arrange(x) |>
  dplyr::group_by(x, .add = TRUE) |>
  dplyr::group_split()

is_even_vec1 <- split_by_unique_x[[1]] |>
  dplyr::arrange(y) |>
  dplyr::mutate(is_even = ifelse(seq_along(y) %% 2 == 0, 1, 0))

is_even_vec2 <- split_by_unique_x[[2]] |>
  dplyr::arrange(y) |>
  dplyr::mutate(is_even = ifelse(seq_along(y) %% 2 == 0, 1, 0))

is_even_vec3 <- split_by_unique_x[[3]] |>
  dplyr::arrange(y) |>
  dplyr::mutate(is_even = ifelse(seq_along(y) %% 2 == 0, 1, 0))

is_even_vec4 <- split_by_unique_x[[4]] |>
  dplyr::arrange(y) |>
  dplyr::mutate(is_even = ifelse(seq_along(y) %% 2 == 0, 1, 0))

is_even_vec5 <- split_by_unique_x[[5]] |>
  dplyr::arrange(y) |>
  dplyr::mutate(is_even = ifelse(seq_along(y) %% 2 == 0, 1, 0))



#box_points$x <- box_points$x + hex_size * ifelse(is_even_vec == 1, 1, 0)
box_points <- bind_rows(is_even_vec1, is_even_vec2, is_even_vec3, is_even_vec4, is_even_vec5)
x_bounds <- box_points$x |> unique()
x_shift <- x_bounds[-length(x_bounds)] - x_bounds[-1]

box_points$x <- box_points$x + x_shift/2 * ifelse(is_even_vec == 1, 0, 1)

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


#box_points <- bind_rows(box_points, middle_box_points)
## Generate all coordinates of hexagons
hex_grid_n <- full_hex_grid(box_points)

ggplot(data = hex_grid_n, aes(x = x, y = y)) + geom_polygon(fill = "white", color = "black", aes(group = id)) +
  geom_point(data = box_points, aes(x = x, y = y), color = "red") +
  geom_point(data = s_curve_noise_umap, aes(x = UMAP1, y = UMAP2), color = "blue")
