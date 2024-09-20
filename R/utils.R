calc_iss <- function(patterns, pattern_semantics) {
  patterns |>
    mutate(
      iss = map_dbl(
        pattern,
        \(pattern) {
          cor(atanh(pattern), pattern_semantics, use = "pairwise")
        }
      ),
      .keep = "unused"
    )
}

calc_iss_stats <- function(data, ..., .by = c(cca_id, time_id)) {
  data |>
    summarise(
      broom::tidy(t.test(iss, ...)),
      .by = {{ .by }}
    )
}

calc_igs <- function(patterns, patterns_group) {
  patterns |>
    mutate(
      igs = map2_dbl(
        pattern, cca_id,
        \(pat, cca) {
          with(
            patterns_group,
            cor(
              pattern[[which(cca_id == cca)]],
              pat,
              use = "pairwise"
            )
          )
        }
      )
    ) |>
    select(!pattern)
}

calc_sync_smc <- function(sync_whole_trials, smc) {
  sync_whole_trials |>
    mutate(
      mantel = map(neu_sync, ~ vegan::mantel(.x, smc)),
      .keep = "unused"
    )
}

calc_clusters_stats <- function(stats, stats_permuted,
                                by = "cca_id",
                                col_statistic = statistic,
                                col_p_value = p.value,
                                col_time_id = time_id,
                                col_id_permuted = starts_with("tar")) {
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
    )
  clusters |>
    left_join(
      clusters_permuted |>
        select(cca_id, cluster_mass_perm = cluster_mass) |>
        chop(cluster_mass_perm),
      by = by
    ) |>
    mutate(
      p_perm = map2_dbl(
        cluster_mass_perm,
        cluster_mass,
        ~ (sum(.x >= .y) + 1) / (length(.x) + 1)
      )
    )
}

find_cluster <- function(statistic, p.value,
                         index = NULL,
                         keep = c("all", "largest"),
                         alpha = 0.05) {
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
        start, end,
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

convert_p2_p1 <- function(statistic, p.value,
                          alternative = c("greater", "less")) {
  alternative <- match.arg(alternative)
  ifelse(
    xor(alternative == "greater", statistic > 0),
    1 - p.value / 2,
    p.value / 2
  )
}

tidy_mantel <- function(mantel) {
  tibble(
    statistic = mantel$statistic,
    p.value = mantel$signif,
    method = mantel$method
  )
}

get_resid <- function(y, x) {
  resid(lm(y ~ x, na.action = na.exclude))
}