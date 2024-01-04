# statistics for subsequent memory effect
extract_stats_sme <- function(dat) {
  dat |>
    group_by(region_id) |>
    group_modify(
      ~ .x |>
        rcompanion::pairwisePermutationMatrix(
          mean_fisher_z ~ response_type_shared,
          data = _,
          method = "holm"
        ) |>
        tidy_pairwise()
    ) |>
    ungroup() |>
    mutate(across(c(x, y), ~ factor(.x, c("Rem", "Know", "Unsure", "New"))))
}

tidy_pairwise <- function(m, name_p_col = "p.value") {
  m[c("Unadjusted", "Adjusted")] |>
    purrr::map(t) |> # the matrix for unadjusted is upper triangle
    purrr::map(stretch, name_value = name_p_col) |>
    bind_rows(.id = "type") |>
    mutate(
      type = case_match(
        type,
        "Unadjusted" ~ "",
        "Adjusted" ~ ".adj"
      )
    ) |>
    pivot_wider(
      names_from = "type",
      names_prefix = name_p_col,
      values_from = all_of(name_p_col)
    )
}

# Extract statistics for each type of prediction
extract_stats_pred_perf <- function(dat, mem_perf) {
  dat |>
    left_join(mem_perf, by = "subj_id", relationship = "many-to-many") |>
    summarise(
      cor.test(mean_fisher_z, score, alternative = "greater") |>
        broom::tidy(),
      .by = c(region_id, window_id, index_name)
    )
}

extract_stats_pred_content <- function(dat, simil_content, ...,
                                       covariate = NULL,
                                       col_rs = mean_fisher_z,
                                       keep_perms = FALSE) {
  dat |>
    mutate(
      map(
        {{ col_rs }},
        ~ stats_mantel(
          as_dist_vec(.x), simil_content, covariate,
          keep_perms = keep_perms
        )
      ) |>
        list_rbind(),
      .keep = "unused"
    )
}

stats_mantel <- function(x, y, z = NULL, ..., keep_perms = FALSE) {
  stats <- if (is.null(z)) {
    vegan::mantel(x, y, ...)
  } else {
    vegan::mantel.partial(x, y, z, ...)
  }
  stats_tbl <- as_tibble_row(stats[c("statistic", "signif")])
  if (keep_perms) {
    stats_tbl$perm <- list(stats$perm)
  }
  stats_tbl
}

# extract cluster-based permutation p value
extract_cluster_p <- function(stats_real, stats_perm, ...,
                              cols_region = region_id,
                              cols_group = NULL,
                              cols_perm = starts_with("tar")) {
  null_distribution <- extract_stats_cluster(
    stats_perm,
    ...,
    by = c({{ cols_region }}, {{cols_group}}, {{ cols_perm }}),
    keep = "largest"
  ) |>
    summarise(
      cluster_mass_perm = max(cluster_mass),
      .by = c({{cols_group}}, {{ cols_perm }})
    ) |>
    select(!{{ cols_perm }}) |>
    chop(cluster_mass_perm)
  cluster_real <- extract_stats_cluster(
    stats_real,
    ...,
    by = c({{ cols_region }}, {{ cols_group }})
  )
  data <- if (is.null(substitute(cols_group))) {
    bind_cols(cluster_real, null_distribution)
  } else {
    cols_group_chr <- tidyselect::eval_select(
      substitute(cols_group),
      cluster_real
    )
    inner_join(cluster_real, null_distribution, by = names(cols_group_chr))
  }
  data |>
    mutate(
      p_perm = map2_dbl(
        cluster_mass, cluster_mass_perm,
        \(real, perm) mean(perm > real)
      ),
      .keep = "unused"
    )
}

extract_cluster_p_rps <- function(file, ...) {
  dat <- R.matlab::readMat(file)
  stats_real <- expand_grid(
    window_id = 1:47,
    region_id = 1:6
  ) |>
    add_column(
      p.value = as.vector(dat$P.real),
      statistic = as.vector(dat$R.real),
      ...
    )
  clusters_p <- stats_real |>
    extract_stats_cluster(by = region_id) |>
    mutate(
      p_perm = map_dbl(
        cluster_mass,
        ~ mean(dat$rsum.max.mean > .x)
      ),
      ...
    )
  lst(stats_real, clusters_p)
}

extract_stats_cluster <- function(stats, by,
                                  col_p_value = p.value,
                                  col_statistic = statistic,
                                  col_window = window_id,
                                  keep = c("all", "largest"),
                                  alpha = 0.05) {
  keep <- match.arg(keep)
  stats |>
    # order is essential for cluster detection
    arrange({{ col_window }}) |>
    reframe(
      find_cluster({{ col_statistic }}, {{ col_p_value }} < alpha, keep),
      .by = {{ by }}
    )
}

find_cluster <- function(statistic, signif, keep) {
  # https://stackoverflow.com/a/43875717/5996475
  rle_signif <- rle(signif)
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
  if (keep == "largest") {
    clusters <- slice_max(clusters, cluster_mass)
  }
  clusters
}
