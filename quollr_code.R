calculate_effective_x_bins <- function(.data, x = "UMAP1", hex_size = NA){

  if (anyNA(.data[[rlang::as_string(rlang::ensym(x))]])) {
    stop("NAs present")
  }

  if (any(is.infinite(.data |> dplyr::pull({{ x }})))) {
    stop("Inf present")
  }



  if (is.na(hex_size)) {
    cell_area <- 1

    ## To compute the radius of the hexagon outer circle
    hex_size <- sqrt(2 * cell_area / sqrt(3))

  } else {
    if ((hex_size <= 0) || (is.infinite(hex_size))) {
      stop("Invalid hex size value")

    }
  }

  ## To compute the range along x-axis
  xwidth <- diff(range(.data |>
                         dplyr::pull({{ x }})))

  horizontal_spacing <- sqrt(3) * hex_size

  num_bins <- ceiling(xwidth/horizontal_spacing)
  num_bins

}

calculate_effective_y_bins <- function(.data, y = "UMAP2", hex_size = NA){

  if (anyNA(.data[[rlang::as_string(rlang::ensym(y))]])) {
    stop("NAs present")
  }

  if (any(is.infinite(.data |> dplyr::pull({{ y }})))) {
    stop("Inf present")
  }

  if (is.na(hex_size)) {
    cell_area <- 1

    ## To compute the radius of the hexagon outer circle
    hex_size <- sqrt(2 * cell_area / sqrt(3))

  } else {
    if ((hex_size <= 0) || (is.infinite(hex_size))) {
      stop("Invalid hex size value")

    }
  }


  ## To compute the range along x-axis
  ywidth <- diff(range(.data |>
                         dplyr::pull({{ y }})))

  vertical_spacing <- 3 * hex_size/2

  num_bins <- ceiling(ywidth/vertical_spacing)
  num_bins


}


# Function to generate even y-values for a given x-value
generate_even_y <- function(data) {
  data <- data |>
    dplyr::arrange(y) |>
    dplyr::mutate(is_even = ifelse(seq_along(y) %% 2 == 0, 1, 0))

  return(data)
}

generate_full_grid_centroids <- function(nldr_df, x = "UMAP1", y = "UMAP2",
                                         num_bins_x, num_bins_y, buffer_size = NA, hex_size = NA){

  ## hex size is not provided
  if (is.na(hex_size)) {
    ## To compute the diameter of the hexagon
    cell_area <- 1

    hex_size <- sqrt(2 * cell_area / sqrt(3))
    message(paste0("hex_size set to ", hex_size, "."))

  }

  ## If number of bins along the x-axis is not given
  if (is.na(num_bins_x)) {
    ## compute the number of bins along the x-axis
    num_bins_x <- calculate_effective_x_bins(.data = nldr_df, x = x, hex_size = hex_size)


  }

  ## If number of bins along the y-axis is not given
  if (is.na(num_bins_y)) {
    num_bins_y <- calculate_effective_y_bins(.data = nldr_df, y = y, hex_size = hex_size)

  }




  ## Buffer size is not provided
  if (is.na(buffer_size)) {
    buffer_size <- hex_size/2
    message(paste0("Buffer set to ", buffer_size, "."))

  } else {

    ## Buffer size is exceeds
    if (buffer_size > (hex_size/2)) {
      stop(paste0("Buffer exceeds than ", hex_size/2, ". Need to assign a value less than ", hex_size/2, "."))

    }


  }

  ## Compute hex grid bound values along the x and y axis
  x_bounds <- seq(min(nldr_df[[rlang::as_string(rlang::sym(x))]]) - buffer_size,
                  max(nldr_df[[rlang::as_string(rlang::sym(x))]]) + buffer_size, length.out = num_bins_x)
  y_bounds <- seq(min(nldr_df[[rlang::as_string(rlang::sym(y))]]) - buffer_size,
                  max(nldr_df[[rlang::as_string(rlang::sym(y))]]) + buffer_size, length.out = num_bins_y)

  ## Generate the all the hex box points
  box_points <- expand.grid(x = x_bounds, y = y_bounds)

  # For each x-value, generate even y-values
  box_points <- box_points |>
    dplyr::arrange(x) |>
    dplyr::group_by(x) |>
    dplyr::group_modify(~ generate_even_y(.x)) |>
    tibble::as_tibble()

  ## Shift the x values of the even rows
  if (length(x_bounds) == 1) {

    if (length(y_bounds) == 1) {
      ## If there is only one bin

      box_points <- tibble::tibble(x = mean(nldr_df[[rlang::as_string(rlang::sym(x))]]),
                                   y = mean(nldr_df[[rlang::as_string(rlang::sym(y))]]))

    } else {
      ## If there is only one bin along x-axis
      x_shift <- 0

      box_points <- box_points |>
        dplyr::select(-is_even)

    }


  } else{

    x_shift <- unique(box_points$x)[2] - unique(box_points$x)[1]

    box_points$x <- box_points$x + x_shift/2 * ifelse(box_points$is_even == 1, 1, 0)

    box_points <- box_points |>
      dplyr::select(-is_even)

  }

  box_points

}

extract_coord_of_shifted_hex_grid <- function(nldr_df, x = "UMAP1", y = "UMAP2",
                                         num_bins_x, num_bins_y, shift_x = 0, shift_y = 0, buffer_size = NA, hex_size = NA){

  ## hex size is not provided
  if (is.na(hex_size)) {
    ## To compute the diameter of the hexagon
    cell_area <- 1

    hex_size <- sqrt(2 * cell_area / sqrt(3))
    message(paste0("hex_size set to ", hex_size, "."))

  }

  ## If number of bins along the x-axis is not given
  if (is.na(num_bins_x)) {
    ## compute the number of bins along the x-axis
    num_bins_x <- calculate_effective_x_bins(.data = nldr_df, x = x, hex_size = hex_size)


  }

  ## If number of bins along the y-axis is not given
  if (is.na(num_bins_y)) {
    num_bins_y <- calculate_effective_y_bins(.data = nldr_df, y = y, hex_size = hex_size)

  }




  ## Buffer size is not provided
  if (is.na(buffer_size)) {
    buffer_size <- hex_size/2
    message(paste0("Buffer set to ", buffer_size, "."))

  } else {

    ## Buffer size is exceeds
    if (buffer_size > (hex_size/2)) {
      stop(paste0("Buffer exceeds than ", hex_size/2, ". Need to assign a value less than ", hex_size/2, "."))

    }


  }

  cell_diameter <- hex_size * 2


  ## Shift is not compatible
  if ((abs(shift_x) > (cell_diameter)) | (abs(shift_y) > (cell_diameter))) {
    stop(paste0("Shifted amount is not compatibel. Need to use a value less than or equal ", cell_diameter, "."))
  }

  x_bounds <- seq(min(nldr_df[[rlang::as_string(rlang::sym(x))]]) - buffer_size,
                  max(nldr_df[[rlang::as_string(rlang::sym(x))]]) + buffer_size, length.out = num_bins_x)
  y_bounds <- seq(min(nldr_df[[rlang::as_string(rlang::sym(y))]]) - buffer_size,
                  max(nldr_df[[rlang::as_string(rlang::sym(y))]]) + buffer_size, length.out = num_bins_y)

  box_points <- expand.grid(x = x_bounds, y = y_bounds)

  # For each x-value, generate even y-values
  box_points <- box_points |>
    dplyr::arrange(x) |>
    dplyr::group_by(x) |>
    dplyr::group_modify(~ generate_even_y(.x)) |>
    tibble::as_tibble()

  ## Shift the x values of the even rows
  if (length(x_bounds) == 1) {
    ## If there is only one bin
    x_shift <- 0

  } else {

    x_shift <- unique(box_points$x)[2] - unique(box_points$x)[1]

  }


  box_points$x <- box_points$x + x_shift/2 * ifelse(box_points$is_even == 1, 1, 0)

  box_points <- box_points |>
    dplyr::select(-is_even)

  box_points

}



extract_hexbin_centroids <- function(hex_full_count_df) {

  df_bin_centroids <- hex_full_count_df[complete.cases(hex_full_count_df[["std_counts"]]), ] |>
    dplyr::select("c_x", "c_y", "hexID", "std_counts") |>
    dplyr::distinct() |>
    dplyr::rename(c("x" = "c_x", "y" = "c_y"))

  return(df_bin_centroids)
}

extract_hexbin_mean <- function(nldr_df_with_hex_id) {

  ## To compute hexagonal bin means
  df_cell_data <- nldr_df_with_hex_id |>
    dplyr::select(-ID) |>
    dplyr::group_by(hb_id) |>
    dplyr::summarise(dplyr::across(tidyselect::everything(), mean))

  ## To rename the columns
  names(df_cell_data) <- c("hexID", "x", "y")

  ## To compute standardise counts
  df_with_std_counts <- compute_std_counts(nldr_df = nldr_df_with_hex_id)

  ## To create the hexbin means info dataset
  df_bin_centroids <- dplyr::inner_join(df_cell_data, df_with_std_counts, by = c("hexID" = "hb_id")) |>
    dplyr::select(x, y, hexID, std_counts)

  return(df_bin_centroids)
}

triangulate_bin_centroids <- function(.data, x = "x", y = "y"){
  tr1 <- tripack::tri.mesh(.data[[x]], .data[[y]])
  return(tr1)
}


generate_edge_info <- function(triangular_object) {

  # Create a data frame with x and y coordinate values from the triangular object
  tr_df <- tibble::tibble(x = triangular_object$x, y = triangular_object$y,
                          ID = 1:length(triangular_object$x))  # Add ID numbers for joining with from and to points in tr_arcs

  # Extract the triangles from the triangular object
  trang <- tripack::triangles(triangular_object)
  trang <- tibble::as_tibble(trang)

  # Create data frames with from-to edges
  tr_arcs_df <- tibble::tibble(from = c(trang$node1, trang$node1, trang$node2),
                               to = c(trang$node2, trang$node3, trang$node3))

  ## To extract unique combinations
  tr_arcs_df <- tr_arcs_df |>
    dplyr::mutate(x = pmin(from, to), y = pmax(from, to)) |>
    dplyr::distinct(x, y) |>
    dplyr::rename(c("from" = "x", "to" = "y"))

  ## Extract from and to values a vectors
  from_vec <- tr_arcs_df$from
  to_vec <- tr_arcs_df$to

  ## Map from and to coordinates
  tr_from_to_df_coord <- dplyr::left_join(tr_arcs_df, tr_df, by = c("from" = "ID")) |>
    dplyr::rename(c("x_from" = "x", "y_from" = "y"))
  tr_from_to_df_coord <- dplyr::left_join(tr_from_to_df_coord, tr_df, by = c("to" = "ID"))|>
    dplyr::rename(c("x_to" = "x", "y_to" = "y"))

  return(tr_from_to_df_coord)


}

cal_2d_dist <- function(tr_from_to_df_coord, start_x = "x_from", start_y = "y_from", end_x = "x_to",
                        end_y = "y_to", select_col_vec = c("from", "to", "distance")) {
  # Calculate the 2D distances
  tr_from_to_df_coord$distance <- lapply(seq(nrow(tr_from_to_df_coord)), function(x) {
    start <- unlist(tr_from_to_df_coord[x, c(start_x, start_y)])
    end <- unlist(tr_from_to_df_coord[x, c(end_x, end_y)])
    sqrt(sum((start - end)^2))
  })

  # Create a data frame with the from-to relationships and distances
  tr_from_to_df_coord <- tr_from_to_df_coord |>
    dplyr::select(tidyselect::all_of(select_col_vec))

  # Convert the distances to a vector and return the data frame
  tr_from_to_df_coord$distance <- unlist(tr_from_to_df_coord$distance, use.names = FALSE)
  return(tr_from_to_df_coord)
}

colour_long_edges <- function(distance_edges, benchmark_value, tr_from_to_df_coord, distance_col) {

  # Create the tibble with x and y coordinates
  tr_df <- tibble::tibble(x = c(tr_from_to_df[["x_from"]], tr_from_to_df[["x_to"]]),
                          y = c(tr_from_to_df[["y_from"]], tr_from_to_df[["y_to"]])) |>
    dplyr::distinct()

  # label small and long edges
  distance_edges <- distance_edges |>
    dplyr::mutate(type = dplyr::if_else(!!as.name(distance_col) < benchmark_value, "small_edges", "long_edges"))

  # Merge edge information with distance data
  tr_from_to_df_coord <- dplyr::inner_join(tr_from_to_df_coord, distance_edges, by = c("from", "to"))

  # Create the triangular mesh plot with colored long edges
  tri_mesh_plot <- ggplot2::ggplot(tr_df, aes(x = x, y = y)) +
    ggplot2::geom_segment(
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to, color = type),
      data = tr_from_to_df_coord
    ) +
    ggplot2::geom_point(size = 1, colour = "#33a02c") +
    ggplot2::coord_equal() +
    ggplot2::scale_colour_manual(values = c("#de2d26", "#636363"))

  return(tri_mesh_plot)
}

remove_long_edges <- function(distance_edges, benchmark_value, tr_from_to_df_coord,
                              distance_col) {
  # Create the tibble with x and y coordinates
  tr_df <- tibble::tibble(x = c(tr_from_to_df[["x_from"]], tr_from_to_df[["x_to"]]),
                          y = c(tr_from_to_df[["y_from"]], tr_from_to_df[["y_to"]])) |>
    dplyr::distinct()

  # Filter small edges
  distance_df_small_edges <- distance_edges |>
    dplyr::filter(!!as.name(distance_col) < benchmark_value)

  # Merge edge information with distance data
  tr_from_to_df_coord <- dplyr::inner_join(tr_from_to_df_coord, distance_df_small_edges,
                                          by = c("from", "to"))


  ## Create the triangular mesh plot after removing the long edges
  tri_mesh_plot <- ggplot2::ggplot(tr_df, aes(x = x, y = y)) + ggplot2::geom_segment(aes(x = x_from,
                                                                                         y = y_from, xend = x_to, yend = y_to), data = tr_from_to_df_coord) +
    ggplot2::geom_point(size = 1, colour = "#33a02c") + ggplot2::coord_equal() + ggplot2::labs(color=NULL)
  return(tri_mesh_plot)

}

# generate_full_grid_centroids <- function(hexdf_data){
#
#   ## Generate initial grid
#   full_centroids1 <- tibble::as_tibble(expand.grid(x = seq(min(hexdf_data$x),max(hexdf_data$x), ggplot2::resolution(hexdf_data$x, FALSE) * 2), y = seq(min(hexdf_data$y),max(hexdf_data$y), ggplot2::resolution(hexdf_data$y, FALSE) * 2)))
#
#   ## Generate shifted grid
#   full_centroids2 <- tibble::tibble(x = full_centroids1$x + ggplot2::resolution(hexdf_data$x, FALSE), y = full_centroids1$y + ggplot2::resolution(hexdf_data$y, FALSE))
#
#   ## Combine all
#   full_centroids <- dplyr::bind_rows(full_centroids1, full_centroids2)
#
#   return(full_centroids)
#
#
# }
## Generate hexagonal coordinates by passing all the centroids
gen_hex_coordinates <- function(all_centroids_df, hex_size = NA){

  # angle_rad_vec <- c()
  #
  # for (i in 1:6) {
  #   angle_deg <- 60 * i - 30
  #   angle_rad = pi / 180 * angle_deg
  #
  #   angle_rad_vec <- append(angle_rad_vec, angle_rad)
  # }

  ## hex size is not provided
  if (is.na(hex_size)) {
    ## To compute the diameter of the hexagon
    cell_area <- 1

    hex_size <- sqrt(2 * cell_area / sqrt(3))
    message(paste0("hex_size set to ", hex_size, "."))

  }

  #angle_rad_vec <- c(0.5235988, 1.5707963, 2.6179939, 3.6651914, 4.7123890, 5.7595865)

  # x_add_factor <- sqrt(3)/2 * (hex_size + hex_size/2)  * cos(angle_rad_vec)
  # y_add_factor <- 3/4 * (hex_size + hex_size/2) * sin(angle_rad_vec)

  # min_value_x <- min(s_curve_noise_umap$UMAP1)
  # max_value_x <- max(s_curve_noise_umap$UMAP1)
  # min_value_y <- min(s_curve_noise_umap$UMAP2)
  # max_value_y <- max(s_curve_noise_umap$UMAP2)

  if (NROW(all_centroids_df) == 1) {

    dx <- hex_size
    dy <- hex_size/ sqrt(3) / 2 * 1.15

  } else {

    dx <- (all_centroids_df$x[2] - all_centroids_df$x[1])
    dy <- (all_centroids_df$y[2] - all_centroids_df$y[1])/ sqrt(3) / 2 * 1.15


  }

  x_add_factor <- c(dx, dx, 0, -dx, -dx, 0)
  y_add_factor <- c(dy, -dy, -2 * dy, -dy, dy, 2 * dy)

  # min_value_x <- min(s_curve_noise_umap$UMAP1)
  # max_value_x <- max(s_curve_noise_umap$UMAP1)
  # min_value_y <- min(s_curve_noise_umap$UMAP2)
  # max_value_y <- max(s_curve_noise_umap$UMAP2)
  #
  # buffer_size <- hex_size/2
  #
  # # Adjust the range of values along the x and y axes with the buffer size
  # x_range_adjusted <- c(min_value_x - buffer_size, max_value_x + buffer_size)
  # y_range_adjusted <- c(min_value_y - buffer_size, max_value_y + buffer_size)
  #
  # # Calculate the spacing between centroids in each direction
  # x_spacing <- (x_range_adjusted[2] - x_range_adjusted[1]) / num_bins_x
  # y_spacing <- (y_range_adjusted[2] - y_range_adjusted[1]) / num_bins_y
  #
  # x_add_factor <- x_spacing  * cos(angle_rad_vec)
  # y_add_factor <- y_spacing * sin(angle_rad_vec)

  all_centroids_df_rep <- all_centroids_df |>
    dplyr::mutate(id = 1:NROW(all_centroids_df)) |>
    dplyr::slice(rep(1:n(), each = 6))

  all_centroids_df_rep_split <- all_centroids_df_rep |>
    dplyr::group_by(id) |>
    dplyr::group_split()

  hex_grid <- data.frame(matrix(ncol = 0, nrow = 0))

  for (id_split in all_centroids_df_rep_split) {

    hex_coord_spec <- tibble::tibble(x = id_split$x + x_add_factor ,
                                     y = id_split$y + y_add_factor,
                                     id = id_split$id)

    hex_grid <- bind_rows(hex_grid, hex_coord_spec)


  }


  return(hex_grid)


}

map_hexbin_id <- function(full_centroid_df) {

  # Create an empty tibble to store all centroids with their coordinates
  full_grid_with_hexbin_id <- tibble::tibble(x = numeric(), y = numeric())

  # Iterate over unique y values
  unique_y <- sort(unique(full_centroid_df$y))

  for(y_val in unique_y){

    ## Filter the data set with specific y value
    specific_y_val_df <- full_centroid_df |>
      dplyr::filter(y == y_val) |>
      dplyr::arrange(x) ## ordered the x values

    full_grid_with_hexbin_id <- dplyr::bind_rows(full_grid_with_hexbin_id, specific_y_val_df)

  }

  ## Add the column with hexagonal bin ID
  full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
    dplyr::mutate(hexID = dplyr::row_number())

  ## Rename columns
  full_grid_with_hexbin_id <- full_grid_with_hexbin_id |>
    dplyr::rename("c_x" = "x",
                  "c_y" = "y")

  return(full_grid_with_hexbin_id)


}

map_polygon_id <- function(full_grid_with_hexbin_id, hex_grid) {

  ## Define a dataset to store polygon id
  full_grid_with_polygon_id <- data.frame(matrix(ncol = 0, nrow = 0))

  unique_hex_id <- unique(full_grid_with_hexbin_id$hexID)
  unique_poly_id <- unique(hex_grid$id)

  for (hb_id_spec in unique_hex_id) {

    ## Filter specific hexagon
    full_grid_with_hexbin_id_filtered <- full_grid_with_hexbin_id |>
      filter(hexID == hb_id_spec)

    for (poly_id_spec in unique_poly_id) {

      ## Filter specific polygon
      hex_grid_filtered <- hex_grid |>
        filter(id == poly_id_spec)

      # Compute the range for x and y coordinates of hex grid
      x_range <- range(hex_grid_filtered$x)
      y_range <- range(hex_grid_filtered$y)

      # Check if each coordinate falls within the range
      status_in_x_range <- full_grid_with_hexbin_id_filtered$c_x >= x_range[1] &
        full_grid_with_hexbin_id_filtered$c_x <= x_range[2]
      status_in_y_range <- full_grid_with_hexbin_id_filtered$c_y >= y_range[1] &
        full_grid_with_hexbin_id_filtered$c_y <= y_range[2]

      # Filter rows where both x and y coordinates are within the range
      filtered_rows <- full_grid_with_hexbin_id_filtered |>
        filter(status_in_x_range & status_in_y_range)

      # Add polygon ID to the filtered rows
      filtered_rows <- mutate(filtered_rows, polygon_id = poly_id_spec)

      # Append filtered rows to the dataset
      full_grid_with_polygon_id <- bind_rows(full_grid_with_polygon_id, filtered_rows)
    }
  }

  return(full_grid_with_polygon_id)
}


assign_data <- function(nldr_df, full_grid_with_hexbin_id) {

  ## Compute distances between nldr coordinates and hex bin centroids
  dist_df <- proxy::dist(as.matrix(nldr_df |>
                                     dplyr::select(-ID)), as.matrix(full_grid_with_hexbin_id |>
                                                                      dplyr::select(-hexID)), method = "Euclidean")

  ## Columns that gives minimum distances
  min_column <- apply(dist_df, 1, which.min)

  nldr_df <- nldr_df |>
    dplyr::mutate(hb_id = full_grid_with_hexbin_id$hexID[min_column])

  return(nldr_df)

}

compute_std_counts <- function(nldr_df) {

  df_with_std_counts <- nldr_df |>
    dplyr::count(hb_id) |>
    dplyr::mutate(std_counts = n/max(n, na.rm = TRUE)) |>
    dplyr::select(-n)

  return(df_with_std_counts)

}

generate_full_grid_info <- function(full_grid_with_polygon_id, df_with_std_counts, hex_grid) {

  ## To assign standardize counts for hex bins
  df_bin_centroids_all <- dplyr::left_join(full_grid_with_polygon_id, df_with_std_counts, by = c("hexID" = "hb_id"))

  ## Since hexagon has 6 coordinates, need to repeat 6 times
  df_bin_centroids_all_rep <- df_bin_centroids_all |>
    dplyr::slice(rep(1:dplyr::n(), each = 6)) |>
    dplyr::arrange(polygon_id)

  ## Join with hexagonal coordinates
  hex_full_count_df <- dplyr::bind_cols(hex_grid, df_bin_centroids_all_rep)

  return(hex_full_count_df)


}


find_pts_in_hexbins <- function(full_grid_with_hexbin_id, nldr_data_with_hb_id) {

  ## Dataframe to store points info
  pts_df <- data.frame(matrix(ncol = 0, nrow = 0))

  for (i in 1:length(nldr_data_with_hb_id$hb_id)) {

    ## Filter a hexagon and find the point within that hexagon
    pts_vec <- nldr_data_with_hb_id |>
      dplyr::filter(hb_id == nldr_data_with_hb_id$hb_id[i]) |>
      dplyr::pull(ID)

    ## Store the hexagon ID with the respective points
    hb_pts <- tibble::tibble(hexID = nldr_data_with_hb_id$hb_id[i], pts = list(pts_vec))

    pts_df <- dplyr::bind_rows(pts_df, hb_pts)

  }

  return(pts_df)

}

find_non_empty_bins <- function(nldr_df_with_id, x = "UMAP1", y = "UMAP2",
                                non_empty_bins, is_bin_centroid = TRUE, hex_size = NA,
                                buffer_size = NA) {

  max_bins_along_axis <- ceiling(sqrt(NROW(nldr_df_with_id)))

  ## Since having 1 bin along x or y-axis is not obvious therefore started from 2
  num_bins_comb_df <- expand.grid(num_bins_x_vec = 2:max_bins_along_axis,
                                  num_bins_y_vec = 2:max_bins_along_axis) |>
    dplyr::mutate(num_x = pmin(num_bins_x_vec, num_bins_y_vec), num_y = pmax(num_bins_x_vec, num_bins_y_vec)) |>
    dplyr::distinct(num_x, num_y)

  num_bins_x <- num_bins_comb_df$num_x[1]
  num_bins_y <- num_bins_comb_df$num_y[1]

  ## hex size is not provided
  if (is.na(hex_size)) {
    ## To compute the diameter of the hexagon
    cell_area <- 1

    hex_size <- sqrt(2 * cell_area / sqrt(3))
    message(paste0("hex_size set to ", hex_size, "."))

  }


  ## Buffer size is not provided
  if (is.na(buffer_size)) {
    buffer_size <- hex_size/2
    message(paste0("Buffer set to ", buffer_size, "."))

  } else {

    ## Buffer size is exceeds
    if (buffer_size > (hex_size/2)) {
      stop(paste0("Buffer exceeds than ", hex_size/2, ". Need to assign a value less than ", hex_size/2, "."))

    }


  }

  ### Generate the full grid

  all_centroids_df <- generate_full_grid_centroids(nldr_df = nldr_df_with_id,
                                                   x = x, y = y,
                                                   num_bins_x = num_bins_x,
                                                   num_bins_y = num_bins_y,
                                                   buffer_size = buffer_size,
                                                   hex_size = hex_size)


  hex_grid <- gen_hex_coordinates(all_centroids_df,
                                  hex_size = hex_size)

  full_grid_with_hexbin_id <- map_hexbin_id(all_centroids_df)

  full_grid_with_polygon_id <- map_polygon_id(full_grid_with_hexbin_id, hex_grid)

  nldr_df_with_hb_id <- assign_data(nldr_df = nldr_df_with_id, full_grid_with_hexbin_id)

  df_with_std_counts <- compute_std_counts(nldr_df = nldr_df_with_hb_id)

  hex_full_count_df <- generate_full_grid_info(full_grid_with_polygon_id,
                                               df_with_std_counts, hex_grid)

  ## Do you need to use bin centroids or bin means?
  if (isTRUE(is_bin_centroid)) {
    ## For bin centroids
    df_bin_centroids <- extract_hexbin_centroids(hex_full_count_df = hex_full_count_df)

  } else {
    ## For bin means
    df_bin_centroids <- extract_hexbin_mean(nldr_df_with_hex_id = nldr_df_with_hb_id)

  }

  num_of_non_empty_bins <- NROW(df_bin_centroids)

  i <- 1

  while (num_of_non_empty_bins < non_empty_bins) {
    i <- i + 1

    num_bins_x <- num_bins_comb_df$num_x[i]
    num_bins_y <- num_bins_comb_df$num_y[i]

    ### Generate the full grid

    all_centroids_df <- generate_full_grid_centroids(nldr_df = nldr_df_with_id,
                                                     x = x, y = y,
                                                     num_bins_x = num_bins_x,
                                                     num_bins_y = num_bins_y,
                                                     buffer_size = buffer_size,
                                                     hex_size = hex_size)


    hex_grid <- gen_hex_coordinates(all_centroids_df,
                                    hex_size = hex_size)

    full_grid_with_hexbin_id <- map_hexbin_id(all_centroids_df)

    full_grid_with_polygon_id <- map_polygon_id(full_grid_with_hexbin_id, hex_grid)

    nldr_df_with_hb_id <- assign_data(nldr_df = nldr_df_with_id, full_grid_with_hexbin_id)

    df_with_std_counts <- compute_std_counts(nldr_df = nldr_df_with_hb_id)

    hex_full_count_df <- generate_full_grid_info(full_grid_with_polygon_id,
                                                 df_with_std_counts, hex_grid)


    ## Do you need to use bin centroids or bin means?
    if (isTRUE(is_bin_centroid)) {
      ## For bin centroids
      df_bin_centroids <- extract_hexbin_centroids(hex_full_count_df = hex_full_count_df)

    } else {
      ## For bin means
      df_bin_centroids <- extract_hexbin_mean(nldr_df_with_hex_id = nldr_df_with_hb_id)

    }

    num_of_non_empty_bins <- NROW(df_bin_centroids)

    if (num_of_non_empty_bins >= non_empty_bins) {
      return(list(num_bins_x = num_bins_x, num_bins_y = num_bins_y))
      break
    } else {
      next
    }
  }
}

avg_highD_data <- function(.data, column_start_text = "x") {
  df_b <- .data |>
    dplyr::select(rsample::starts_with(column_start_text), hb_id) |>
    dplyr::group_by(hb_id) |>
    dplyr::summarise(dplyr::across(tidyselect::everything(), mean))

  return(df_b)
}

compute_weights <- function(nldr_df_with_hb_id) {

  ## To get the average of each bin
  bin_val_hexagons <- nldr_df_with_hb_id |>
    dplyr::select(-ID) |>
    dplyr::group_by(hb_id) |>
    dplyr::summarise(dplyr::across(tidyselect::everything(), mean))

  new_col <- paste0("avg_", names(nldr_df_with_hb_id)[1:2] |> tolower())

  names(bin_val_hexagons) <- append("hb_id", new_col)

  ## To calculate distances from average point

  nldr_with_avg_all <- dplyr::inner_join(bin_val_hexagons , nldr_df_with_hb_id,
                                         by = c("hb_id" = "hb_id")) |>
    dplyr::select(-ID)


  nldr_with_avg_all_split <- nldr_with_avg_all |>
    dplyr::group_by(hb_id) |>
    dplyr::group_split()

  col_names1 <- append(names(bin_val_hexagons), names(nldr_df_with_hb_id)[1:2])
  col_names <- append(col_names1, "distance")

  vec <- stats::setNames(1:6, col_names)
  weight_df <- dplyr::bind_rows(vec)[0, ]

  for(i in 1:length(nldr_with_avg_all_split)){

    weighted_mean_df <- nldr_with_avg_all_split[[i]] |> ## These are the weights for weighted mean
      cal_2d_dist(start_x = new_col[1], start_y = new_col[2], end_x = names(nldr_df_with_hb_id)[1],
                  end_y = names(nldr_df_with_hb_id)[2], select_col_vec = col_names)

    weight_df <- dplyr::bind_rows(weight_df, weighted_mean_df)

  }

  return(weight_df)

}

weighted_high_d_data <- function(training_data, nldr_df_with_hb_id, column_start_text = "x") {

  df_all <- dplyr::bind_cols(training_data |> dplyr::select(-ID), nldr_df_with_hb_id)

  weight_df <- compute_weights(nldr_df_with_hb_id = nldr_df_with_hb_id)

  joined_col_names <- names(nldr_df_with_hb_id)[1:2]

  weighted_mean_all <- dplyr::inner_join(df_all, weight_df, by = c("hb_id" = "hb_id",
                                                                   stats::setNames(joined_col_names, joined_col_names))) |>
    mutate(distance_trans =  1/ (distance + 0.05))

  weighted_mean_df_list <- list()

  for (j in 1:(NCOL(training_data |> dplyr::select(-ID)))) {

    weighted_mean_df_list[[j]] <- weighted_mean_all |>
      dplyr::select(hb_id, names(training_data |> dplyr::select(-ID))[j], distance_trans) |>
      dplyr::group_by(hb_id) |>
      dplyr::summarise(dplyr::across(names(training_data |> dplyr::select(-ID))[j], ~ weighted.mean(., distance_trans)))

  }

  weighted_mean <- Reduce(function(dtf1,dtf2) dplyr::full_join(dtf1,dtf2,by="hb_id"),
                          weighted_mean_df_list)


  ## Column names start with x
  weighted_mean <- weighted_mean |>
    dplyr::select(hb_id, tidyselect::starts_with(column_start_text))

  return(weighted_mean)
}

show_langevitour <- function(df, df_b, df_b_with_center_data, benchmark_value = NA,
                             distance_df, distance_col, use_default_benchmark_val = FALSE, column_start_text = "x") {

  ### Define type column
  df <- df |>
    dplyr::select(tidyselect::starts_with(column_start_text)) |>
    dplyr::mutate(type = "data") ## original dataset

  df_b <- df_b |>
    dplyr::filter(hb_id %in% df_b_with_center_data$hexID) |>
    dplyr::mutate(type = "model") ## Data with summarized mean

  ## Reorder the rows of df_b according to the hexID order in df_b_with_center_data
  df_b <- df_b[match(df_b_with_center_data$hexID, df_b$hb_id),] |>
    dplyr::select(-hb_id)

  df_exe <- dplyr::bind_rows(df_b, df)


  if(is.na(benchmark_value)){

    if (isFALSE(use_default_benchmark_val)) {

      tr1 <- triangulate_bin_centroids(df_b_with_center_data, x = "x", y = "y")
      tr_from_to_df <- generate_edge_info(triangular_object = tr1)

      langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = tr_from_to_df$from,
                               lineTo = tr_from_to_df$to, group = df_exe$type, pointSize = 3,
                               levelColors = c("#6a3d9a", "#33a02c"))

    } else {

      benchmark_value <- find_benchmark_value(distance_edges = distance_df, distance_col = distance_col)

      ## Set the maximum difference as the criteria
      distance_df_small_edges <- distance_df |>
        dplyr::filter(!!as.name(distance_col) < benchmark_value)
      ## Since erase brushing is considerd.

      langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from,
                               lineTo = distance_df_small_edges$to, group = df_exe$type, pointSize = 3,
                               levelColors = c("#6a3d9a", "#33a02c"))

    }

  } else {

    ## Check benchmark value is an accepted one
    if (benchmark_value < min(distance_df[[distance_col]])) {
      stop("Benchmark value to remove long edges is too small.")

    }

    if (benchmark_value > max(distance_df[[distance_col]])) {
      stop("Benchmark value to remove long edges is too large.")

    }

    if (isTRUE(use_default_benchmark_val)) {
      stop("Need to set `benchmark_value = NA`.")
    }

    ## Set the maximum difference as the criteria
    distance_df_small_edges <- distance_df |>
      dplyr::filter((!!as.name(distance_col)) < benchmark_value)
    ## Since erase brushing is considerd.

    langevitour::langevitour(df_exe[1:(length(df_exe)-1)], lineFrom = distance_df_small_edges$from,
                             lineTo = distance_df_small_edges$to, group = df_exe$type, pointSize = 3,
                             levelColors = c("#6a3d9a", "#33a02c"))

  }


}

fit_high_d_model <- function(training_data, nldr_df_with_id, x = "UMAP1",
                             y = "UMAP1", num_bins_x = NA, num_bins_y = NA,
                             hex_size = NA, buffer_size = NA,
                             is_bin_centroid = TRUE,
                             is_rm_lwd_hex = FALSE,
                             benchmark_to_rm_lwd_hex = NA,
                             is_avg_high_d = TRUE, column_start_text = "x") {

  ## hex size is not provided
  if (is.na(hex_size)) {
    ## To compute the diameter of the hexagon
    cell_area <- 1

    hex_size <- sqrt(2 * cell_area / sqrt(3))
    message(paste0("hex_size set to ", hex_size, "."))

  }

  ## If number of bins along the x-axis is not given
  if (is.na(num_bins_x)) {
    ## compute the number of bins along the x-axis
    num_bins_x <- calculate_effective_x_bins(.data = nldr_df_with_id, x = x, hex_size = hex_size)


  }

  ## If number of bins along the y-axis is not given
  if (is.na(num_bins_y)) {
    num_bins_y <- calculate_effective_y_bins(.data = nldr_df_with_id, y = y, hex_size = hex_size)

  }

  ## Buffer size is not provided
  if (is.na(buffer_size)) {
    buffer_size <- hex_size/2
    message(paste0("Buffer set to ", buffer_size, "."))

  } else {

    ## Buffer size is exceeds
    if (buffer_size > (hex_size/2)) {
      stop(paste0("Buffer exceeds than ", hex_size/2, ". Need to assign a value less than ", hex_size/2, "."))

    }


  }

  ### Generate the full grid

  all_centroids_df <- generate_full_grid_centroids(nldr_df = nldr_df_with_id,
                                                   x = x, y = y,
                                                   num_bins_x = num_bins_x,
                                                   num_bins_y = num_bins_y,
                                                   buffer_size = buffer_size, hex_size = hex_size)


  hex_grid <- gen_hex_coordinates(all_centroids_df,
                                  hex_size = hex_size)

  full_grid_with_hexbin_id <- map_hexbin_id(all_centroids_df)

  full_grid_with_polygon_id <- map_polygon_id(full_grid_with_hexbin_id, hex_grid)

  nldr_df_with_hb_id <- assign_data(nldr_df = nldr_df_with_id, full_grid_with_hexbin_id)

  df_with_std_counts <- compute_std_counts(nldr_df = nldr_df_with_hb_id)

  hex_full_count_df <- generate_full_grid_info(full_grid_with_polygon_id, df_with_std_counts, hex_grid)

  ## Do you need to use bin centroids or bin means?
  if (isTRUE(is_bin_centroid)) {
    ## For bin centroids
    df_bin_centroids <- extract_hexbin_centroids(hex_full_count_df = hex_full_count_df)

  } else {
    ## For bin means
    df_bin_centroids <- extract_hexbin_mean(nldr_df_with_hex_id = nldr_df_with_hb_id)

  }



  ## Do you need to remove low density hexagons?
  if (isTRUE(is_rm_lwd_hex)) {

    ## if the benchmark value to remove low density hexagons is not provided
    if (is.na(benchmark_to_rm_lwd_hex)) {
      ## first quartile used as the default
      benchmark_to_rm_lwd_hex <- stats::quantile(df_bin_centroids$std_counts,
                                                 probs = c(0,0.25,0.5,0.75,1), names = FALSE)[2]
    }

    ## To identify low density hexagons
    df_bin_centroids_low <- df_bin_centroids |>
      dplyr::filter(std_counts <= benchmark_to_rm_lwd_hex)

    ## To identify low-density hexagons needed to remove by investigating neighbouring mean density
    identify_rm_bins <- find_low_density_hexagons(df_bin_centroids_all = df_bin_centroids,
                                                  num_bins_x = num_bins_x,
                                                  df_bin_centroids_low = df_bin_centroids_low)

    ## To remove low-density hexagons
    df_bin_centroids <- df_bin_centroids |>
      filter(!(hexID %in% identify_rm_bins))



  }


  ## To generate a data set with high-D and 2D training data
  df_all <- dplyr::bind_cols(training_data |> dplyr::select(-ID), nldr_df_with_hb_id)

  ## Do you need to use bin centroids or bin means?
  if (isTRUE(is_avg_high_d)) {

    ## averaged high-D data
    df_bin <- avg_highD_data(.data = df_all, column_start_text = column_start_text)


  } else {

    ## weighted averaged high-D data
    df_bin <- weighted_high_d_data(training_data = training_data,
                                   nldr_df_with_hb_id = nldr_df_with_hb_id,
                                   column_start_text = column_start_text)

  }

  ## high-D model only contains the bins in 2D
  df_bin <- df_bin |>
    dplyr::filter(hb_id %in% df_bin_centroids$hexID)

  return(list(df_bin = df_bin, df_bin_centroids = df_bin_centroids))

}

find_benchmark_value <- function(distance_edges, distance_col) {

  if (any(is.na(distance_edges[[distance_col]]))) {
    stop("NAs present")
  }

  distance_edges <- distance_edges |>
    dplyr::select(!!rlang::sym(distance_col)) |>
    dplyr::mutate(dplyr::across({
      {
        distance_col
      }
    }, \(x) round(x, 3))) |>
    dplyr::arrange(!!rlang::sym(distance_col)) |>  ## Sort the distances
    dplyr::distinct()  ## Find unique distances

  ## Calculate differences between unique distance

  distance_edges <- distance_edges |>
    dplyr::mutate(difference = append(0, apply(distance_edges, 2, diff))) |>
    dplyr::mutate(dplyr::across(difference, ~ round(., 4)))  ## For simplicity

  benchmark_value_vec <- c()

  ## To find the first largest difference (Define a benchmark value
  ## to remove long edges)
  for (i in 1:dim(distance_edges)[1]) {
    if(!is.na(distance_edges$difference[i + 1])){
      if (distance_edges$difference[i] > distance_edges$difference[i + 1]) {
        if (!(is.na(distance_edges$difference[i]))) {
          benchmark_value_vec[i] <- distance_edges$difference[i]
          break
        }
      }
    }
  }

  benchmark_value <- distance_edges[which(distance_edges$difference == benchmark_value_vec[!(is.na(benchmark_value_vec))]),
                                    1] |>  # To get the first value which contain large difference
    dplyr::pull(distance) |>
    dplyr::nth(1)


  if (is.na(benchmark_value)) {
    ## first quartile used as the default
    benchmark_value <- stats::quantile(distance_edges$distance,
                                               probs = c(0,0.25,0.5,0.75,1), names = FALSE)[2]

  }

  benchmark_value


}

compute_mean_density_hex <- function(df_bin_centroids, num_bins_x) {

  # To store mean densities of hexagons
  mean_density_vec <- c()

  for (i in 1:length(df_bin_centroids$hexID)) {

    ## Identify neighbors of a specific hex bin
    neighbor_df <- df_bin_centroids |>
      dplyr::filter((hexID == (df_bin_centroids$hexID[i] + 1)) | (hexID == (df_bin_centroids$hexID[i] - 1)) |
                      (hexID == (df_bin_centroids$hexID[i] + (num_bins_x + 1))) |
                      (hexID == (df_bin_centroids$hexID[i] + num_bins_x)) |
                      (hexID == (df_bin_centroids$hexID[i] - (num_bins_x + 1))) |
                      (hexID == (df_bin_centroids$hexID[i] - num_bins_x)))

    mean_density <- neighbor_df |>
      dplyr::pull(std_counts) |>
      sum()/NROW(neighbor_df) ## The reason to take the mean is to check the density in a considerable amount

    mean_density_vec <- append(mean_density_vec, mean_density)

  }

  df_bin_centroids <- df_bin_centroids |>
    dplyr::mutate(mean_density = mean_density_vec)

  return(df_bin_centroids)

}

find_low_density_hexagons <- function(df_bin_centroids_all, num_bins_x, df_bin_centroids_low) {
  ## To compute mean density of hexagons
  df_bin_centroids <- compute_mean_density_hex(df_bin_centroids_all, num_bins_x)
  mean_density_vec <- df_bin_centroids$mean_density

  df_bin_centroids_low <- df_bin_centroids |>
    dplyr::filter(hexID %in% df_bin_centroids_low$hexID)

  ## Take first quartile
  benchmark_mean_dens_rm_hex <- stats::quantile(mean_density_vec,
                                                probs = c(0,0.25,0.5,0.75,1),
                                                names = FALSE)[2]

  remove_bins <- c()

  ## Check only already identified low-density hexagons
  for (i in 1:length(df_bin_centroids_low$hexID)) {

    df_bin_centroids_coordinates_spec_bin <- df_bin_centroids_low |>
      dplyr::filter(hexID == df_bin_centroids_low$hexID[i])

    bin_ID <- df_bin_centroids_coordinates_spec_bin |>
      dplyr::pull(hexID)


    if(df_bin_centroids_coordinates_spec_bin$mean_density < benchmark_mean_dens_rm_hex){
      remove_bins <- append(remove_bins, bin_ID)
    }
  }

  return(remove_bins)
}

compute_aic <- function(p, total, num_bins, num_obs) {
  mse <- mean(total)
  aic <- 2*num_bins*p + num_obs*p*log(mse)
  return(aic)
}

generate_summary <- function(test_data, prediction_df, df_bin, col_start = "x") {

  ## Rename columns to avoid conflicts
  names(df_bin)[-1] <- paste0("model_high_d_", names(df_bin)[-1])

  prediction_df <- prediction_df |>
    dplyr::left_join(df_bin, by = c("pred_hb_id" = "hb_id")) ## Map high-D averaged/weighted mean coordinates

  prediction_df <- prediction_df |>
    dplyr::left_join(training_data, by = c("ID" = "ID")) ## Map high-D data

  cols <- paste0(col_start, 1:(NCOL(df_bin) - 1))
  high_d_model_cols <- paste0("model_high_d_", col_start, 1:(NCOL(df_bin) - 1))
  error_cols <- paste0("error_square_", col_start, 1:(NCOL(df_bin) - 1))

  summary_df <- (prediction_df[, cols] - prediction_df[, high_d_model_cols])^2
  names(summary_df) <- error_cols

  row_wise_total_error <- rowSums(summary_df[, error_cols, drop = FALSE])

  aic <-  compute_aic((NCOL(df_bin) - 1), row_wise_total_error,
                    NROW(df_bin), NROW(training_data))
  mse <-  mean(row_wise_total_error)

  return(list(mse = mse, aic = aic))

}

predict_2d_embeddings <- function(test_data, df_bin_centroids, df_bin, type_NLDR = "UMAP") {

  ## Compute distances between nldr coordinates and hex bin centroids
  dist_df <- proxy::dist(as.matrix(test_data |> dplyr::select(-ID)),
                         as.matrix(df_bin |> dplyr::select(-hb_id)), method = "Euclidean")

  ## Columns that gives minimum distances
  min_column <- apply(dist_df, 1, which.min)

  test_data <- test_data |>
    dplyr::mutate(pred_hb_id = df_bin$hb_id[min_column])

  ## Obtain 2D coordinate of the nearest high-D centroid
  match_indices <- match(test_data$pred_hb_id, df_bin_centroids$hexID)
  predict_centroid_coord_2D <- dplyr::bind_cols(test_data, emb_1 = df_bin_centroids$x[match_indices],
                                                emb_2 = df_bin_centroids$y[match_indices])

  predict_centroid_coord_2D <- predict_centroid_coord_2D |>
    dplyr::select(emb_1, emb_2, ID, pred_hb_id)

  ## Rename columns
  names(predict_centroid_coord_2D) <- c(paste0("pred_", type_NLDR, "_", 1:2), "ID", "pred_hb_id")

  return(predict_centroid_coord_2D)

}

geom_trimesh <- function(mapping = NULL, data = NULL, stat = "trimesh",
                         position = "identity", show.legend = NA, na.rm = FALSE, inherit.aes = TRUE,
                         ...) {
  ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = GeomTrimesh,
                 position = position, show.legend = show.legend, inherit.aes = inherit.aes,
                 params = list(na.rm = na.rm, ...))
}

GeomTrimesh <- ggplot2::ggproto("GeomTrimesh",
                                ggplot2::Geom,
                                required_aes = c("x", "y", "xend", "yend"),
                                default_aes = ggplot2::aes(
                                  shape = 19,
                                  linetype = 1,
                                  linewidth = 0.5,
                                  size = 0.5,
                                  alpha = NA,
                                  colour = "#33a02c"
                                ),
                                draw_key = ggplot2::draw_key_point,
                                draw_panel = function(data, panel_scales, coord) {

                                  vertices <- tibble::tibble(
                                    x = data$x,
                                    y = data$y,
                                    colour = rep("#33a02c", nrow(data)),
                                    shape = data$shape,
                                    size = rep(2, nrow(data)),
                                    fill = rep("#33a02c", nrow(data)),
                                    alpha = data$alpha,
                                    stroke = 0.5,
                                    stringsAsFactors = FALSE
                                  )

                                  trimesh <- tibble::tibble(
                                    x = data$x,
                                    xend = data$xend,
                                    y = data$y,
                                    yend = data$yend,
                                    PANEL = data$PANEL,
                                    group = data$group,
                                    size = data$size,
                                    linetype = data$linetype,
                                    linewidth = data$linewidth,
                                    alpha = data$alpha,
                                    colour = data$colour
                                  )

                                  ggplot2:::ggname(
                                    "geom_trimesh",
                                    grid::grobTree(
                                      ggplot2::GeomSegment$draw_panel(trimesh, panel_scales, coord),
                                      ggplot2::GeomPoint$draw_panel(vertices, panel_scales, coord)
                                    )
                                  )
                                }
)

stat_trimesh <- function(mapping = NULL, data = NULL, geom = GeomTrimesh$default_aes(),
                         position = "identity", show.legend = NA, outliers = TRUE, inherit.aes = TRUE,
                         ...) {
  ggplot2::layer(
    stat = StatTrimesh,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(outliers = outliers, ...)
  )
}

StatTrimesh <- ggplot2::ggproto(
  "StatTrimesh",
  ggplot2::Stat,
  compute_group = function(data, scales, outliers = TRUE) {
    tr1 <- tripack::tri.mesh(data$x, data$y, duplicate = "remove")
    tr_df <- tibble::tibble(x = tr1$x, y = tr1$y)  ## Create a dataframe with tri.mesh x and y coordinate values
    tr_df <- tr_df |>
      dplyr::mutate(ID = dplyr::row_number())  ## To add ID numbers, because to join with from and to points in tri$arcs

    trang <- tripack::triangles(tr1)
    trang <- tibble::as_tibble(trang)

    tr_arcs_df1 <- tibble::tibble(from = trang$node1, to = trang$node2)  ## Create dataframe with from and to edges
    tr_arcs_df2 <- tibble::tibble(from = trang$node1, to = trang$node3)
    tr_arcs_df3 <- tibble::tibble(from = trang$node2, to = trang$node3)
    tr_arcs_df <- dplyr::bind_rows(tr_arcs_df1, tr_arcs_df2, tr_arcs_df3)  ## Create dataframe with from and to edges

    ## To obtain x and values of from to in a dataframe
    vec <- stats::setNames(rep("", 6), c("from", "to", "x_from", "y_from",
                                         "x_to", "y_to"))  ## Define column names
    # Initialize an empty dataframe to store data in a specific
    # format
    tr_from_to_df_coord <- dplyr::bind_rows(vec)[0, ]
    tr_from_to_df_coord <- tr_from_to_df_coord |>
      dplyr::mutate_if(is.character, as.numeric)

    for (i in 1:NROW(tr_arcs_df)) {
      from_row <- tr_df |>
        dplyr::filter(dplyr::row_number() == (tr_arcs_df |>
                                                dplyr::pull(from) |>
                                                dplyr::nth(i)))
      to_row <- tr_df |>
        dplyr::filter(dplyr::row_number() == (tr_arcs_df |>
                                                dplyr::pull(to) |>
                                                dplyr::nth(i)))
      tr_from_to_df_coord <- tr_from_to_df_coord |>
        tibble::add_row(from = from_row |>
                          dplyr::pull(ID),
                        to = to_row |>
                          dplyr::pull(ID),
                        x_from = from_row |>
                          dplyr::pull(x),
                        y_from = from_row |>
                          dplyr::pull(y),
                        x_to = to_row |>
                          dplyr::pull(x),
                        y_to = to_row |>
                          dplyr::pull(y))  ## Add vector as an appending row to the dataframe
    }

    trimesh <- tibble::tibble(x = tr_from_to_df_coord$x_from,
                              y = tr_from_to_df_coord$y_from,
                              xend = tr_from_to_df_coord$x_to,
                              yend = tr_from_to_df_coord$y_to,
                              PANEL = as.factor(rep(1, nrow(tr_from_to_df_coord))),
                              group = rep(-1, nrow(tr_from_to_df_coord)),
                              size = rep(0.5, nrow(tr_from_to_df_coord)),
                              linetype = rep(1, nrow(tr_from_to_df_coord)),
                              linewidth = rep(0.5, nrow(tr_from_to_df_coord)),
                              alpha = rep(NA, nrow(tr_from_to_df_coord)),
                              colour = rep("#636363", nrow(tr_from_to_df_coord)))
    trimesh
  },
  required_aes = c("x", "y")
)
