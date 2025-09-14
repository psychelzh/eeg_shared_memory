calc_clusters_stats <- function(
  stats,
  stats_permuted,
  by = "cca_id",
  col_statistic = statistic,
  col_p_value = p.value,
  col_time_id = time_id,
  col_id_permuted = .rep,
  alternative = c("greater", "less")
) {
  operator <- switch(match.arg(alternative), greater = `>=`, less = `<=`)
  clusters <- stats |>
    reframe(
      find_cluster(
        {{ col_statistic }},
        {{ col_p_value }},
        {{ col_time_id }}
      ),
      .by = all_of(by)
    )
  clusters_permuted <- stats_permuted |>
    reframe(
      find_cluster(
        {{ col_statistic }},
        {{ col_p_value }},
        {{ col_time_id }},
        keep = "largest"
      ),
      .by = c(all_of(by), {{ col_id_permuted }})
    ) |>
    select(all_of(by), cluster_mass_perm = cluster_mass) |>
    chop(cluster_mass_perm)
  clusters_combined <- if (is_empty(by)) {
    cross_join(clusters, clusters_permuted)
  } else {
    inner_join(clusters, clusters_permuted, by = by)
  }
  clusters_combined |>
    mutate(
      p_perm = map2_dbl(
        cluster_mass_perm,
        cluster_mass,
        ~ (sum(operator(.x, .y)) + 1) / (length(.x) + 1)
      )
    )
}

find_cluster <- function(
  statistic,
  p.value,
  index = NULL,
  keep = c("all", "largest"),
  alpha = 0.05
) {
  keep <- match.arg(keep)
  # https://stackoverflow.com/a/43875717/5996475
  rle_signif <- rle(p.value < alpha)
  if (!any(rle_signif$values)) {
    return(tibble(start = NA, end = NA, cluster_mass = 0))
  }
  end <- cumsum(rle_signif$lengths)
  start <- c(1, lag(end)[-1] + 1)
  clusters <- tibble(
    start = start[rle_signif$values],
    end = end[rle_signif$values]
  ) |>
    mutate(
      cluster_mass = map2_dbl(
        start,
        end,
        \(start, end) {
          sum(statistic[start:end])
        }
      )
    )
  if (!is.null(index)) {
    clusters$start <- index[clusters$start]
    clusters$end <- index[clusters$end]
  }
  if (keep == "largest") {
    clusters <- slice_max(clusters, cluster_mass)
  }
  clusters
}
