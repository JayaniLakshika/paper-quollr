## This script is to generate highd and 2D distances

library(tidyverse)
library(proxy)

data <- read_rds("data/limb_muscles/facs_limb_muscles_pcs_10.rds")

## Compute highd distances
dist_vec <- proxy::dist(data[, 1:10], method = "Euclidean") |> as.vector()

from_vec <- c()
to_vec <- c()
num_obs <- 1:(NROW(data) - 1)

for (obs in num_obs) {

  from_val <- rep(obs, (NROW(data) - obs))
  if ((obs + 1) <= NROW(data)) {
    to_val <- (obs + 1):NROW(data)
  }
  from_vec <- append(from_vec, from_val)
  to_vec <- append(to_vec, to_val)

}

dist_highd <- tibble::tibble(from = from_vec, to = to_vec, dist_highd = dist_vec)

#### Compute 2D distances (tsne author)
tsne_limb <- read_rds("data/limb_muscles/facs_limb_muscles_tsne_perplexity_30.rds")

dist_vec <- proxy::dist(x = tsne_limb[, 1:2], method = "Euclidean") |> as.vector()

from_vec <- c()
to_vec <- c()
num_obs <- 1:(NROW(tsne_limb) - 1)

for (obs in num_obs) {

  from_val <- rep(obs, (NROW(tsne_limb) - obs))
  if ((obs + 1) <= NROW(tsne_limb)) {
    to_val <- (obs + 1):NROW(tsne_limb)
  }
  from_vec <- append(from_vec, from_val)
  to_vec <- append(to_vec, to_val)

}

dist_2d1 <- tibble::tibble(from = from_vec, to = to_vec, dist_2d_tsne_30 = dist_vec) |>
  dplyr::select(dist_2d_tsne_30)

#### Compute 2D distances (tsne best)
tsne_limb <- read_rds("data/limb_muscles/facs_limb_muscles_tsne_perplexity_15.rds")

dist_vec <- proxy::dist(x = tsne_limb[, 1:2], method = "Euclidean") |> as.vector()

from_vec <- c()
to_vec <- c()
num_obs <- 1:(NROW(tsne_limb) - 1)

for (obs in num_obs) {

  from_val <- rep(obs, (NROW(tsne_limb) - obs))
  if ((obs + 1) <= NROW(tsne_limb)) {
    to_val <- (obs + 1):NROW(tsne_limb)
  }
  from_vec <- append(from_vec, from_val)
  to_vec <- append(to_vec, to_val)

}

dist_2d2 <- tibble::tibble(from = from_vec, to = to_vec, dist_2d_tsne_15 = dist_vec) |>
  dplyr::select(dist_2d_tsne_15)

dist_df <- dplyr::bind_cols(dist_2d1, dist_2d2, dist_highd)
dist_df <- dist_df |>
  dplyr::select(from, to, dist_2d_tsne_30, dist_2d_tsne_15, dist_highd)

### Run only once
write_rds(dist_df, "data/limb_muscles/facs_limb_dist_df.rds")

