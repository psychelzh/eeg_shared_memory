calc_mantel <- function(
  data,
  ydis,
  col_xdis = last_col(),
  seed = NULL,
  permutations = 9999L
) {
  force(ydis)
  row_seeds <- if (is.null(seed)) {
    rep(list(NULL), nrow(data))
  } else {
    as.list(partition_row_seeds(data, seed))
  }
  data |>
    mutate(
      mantel = map2(
        pick({{ col_xdis }})[[1]],
        row_seeds,
        \(xdis, row_seed) {
          if (is.null(row_seed)) {
            vegan::mantel(
              xdis,
              ydis,
              permutations = permutations,
              parallel = 1L
            )
          } else {
            with_permutation_seed(
              row_seed,
              vegan::mantel(
                xdis,
                ydis,
                permutations = permutations,
                parallel = 1L
              )
            )
          }
        }
      ),
      .keep = "unused"
    )
}

calc_mantel_partial <- function(
  data,
  ydis,
  zdis,
  col_xdis = last_col(),
  seed = NULL,
  permutations = 9999L
) {
  force(ydis)
  force(zdis)
  row_seeds <- if (is.null(seed)) {
    rep(list(NULL), nrow(data))
  } else {
    as.list(partition_row_seeds(data, seed))
  }
  data |>
    mutate(
      mantel = map2(
        pick({{ col_xdis }})[[1]],
        row_seeds,
        \(xdis, row_seed) {
          if (is.null(row_seed)) {
            vegan::mantel.partial(
              xdis,
              ydis,
              zdis,
              permutations = permutations,
              parallel = 1L
            )
          } else {
            with_permutation_seed(
              row_seed,
              vegan::mantel.partial(
                xdis,
                ydis,
                zdis,
                permutations = permutations,
                parallel = 1L
              )
            )
          }
        }
      ),
      .keep = "unused"
    )
}

extract_stats_mantel <- function(data_mantel) {
  data_mantel |>
    mutate(
      map(mantel, tidy_mantel) |>
        list_rbind(),
      .keep = "unused"
    )
}

tidy_mantel <- function(mantel) {
  tibble(
    statistic = mantel$statistic,
    p.value = mantel$signif,
    method = mantel$method
  )
}
