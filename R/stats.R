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

extract_stats_pred_content <- function(dat, simil_content, col_rs) {
  dat |>
    mutate(
      map(
        {{ col_rs }},
        ~ vegan::mantel(as_dist_vec(.x), simil_content) |>
          unclass() |>
          select_list(all_of(c("statistic", "signif"))) |>
          as_tibble_row()
      ) |>
        list_rbind(),
      .keep = "unused"
    )
}

extract_stats_pred_content_partial <- function(dat, simil_content, covariate,
                                               col_rs) {
  dat |>
    mutate(
      map(
        {{ col_rs }},
        ~ vegan::mantel.partial(as_dist_vec(.x), simil_content, covariate) |>
          unclass() |>
          select_list(all_of(c("statistic", "signif"))) |>
          as_tibble_row()
      ) |>
        list_rbind(),
      .keep = "unused"
    )
}


# extract cluster-based permutation p value
extract_cluster_p <- function(stats_cluster,
                              stats_cluster_perm,
                              col_stats = sum_t,
                              col_start = start,
                              col_end = end) {
  cols_group <- tidyselect::eval_select(
    rlang::expr(!c({{ col_stats }}, {{ col_start }}, {{ col_end }})),
    stats_cluster
  )
  stats_cluster_perm |>
    nest(.by = all_of(names(cols_group))) |>
    left_join(stats_cluster, by = names(cols_group)) |>
    mutate(
      p_perm = map2_dbl(
        data, {{ col_stats }},
        \(perms, real) mean(pull(perms, {{ col_stats }}) > real)
      ),
      .keep = "unused"
    )
}

extract_stats_cluster <- function(stats, .by,
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
      .by = {{ .by }}
    )
}

find_cluster <- function(statistic, signif, keep) {
  # https://stackoverflow.com/a/43875717/5996475
  rle_signif <- rle(signif)
  if (!any(rle_signif$values)) {
    return(tibble(start = NA, end = NA, sum_t = 0))
  }
  end <- cumsum(rle_signif$lengths)
  start <- c(1, lag(end)[-1] + 1)
  clusters <- tibble(
    start = start[rle_signif$values],
    end = end[rle_signif$values]
  ) |>
    mutate(
      sum_t = map2_dbl(
        start, end,
        \(start, end) {
          sum(statistic[start:end])
        }
      )
    )
  if (keep == "largest") {
    clusters <- slice_max(clusters, sum_t)
  }
  clusters
}
