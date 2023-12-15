extract_stats_pred_perf <- function(file_parquet, mem_perf) {
  arrow::open_dataset(file_parquet) |>
    # here we must average across different trials
    summarise(
      mean_fisher_z = mean(fisher_z, na.rm = TRUE),
      .by = c(subj_id, region_id, window_id)
    ) |>
    left_join(mem_perf, by = "subj_id", relationship = "many-to-many") |>
    collect() |>
    summarise(
      cor.test(mean_fisher_z, dprime, alternative = "greater") |>
        broom::tidy(),
      .by = c(region_id, window_id, mem_type)
    )
}

# permutate subject id to get surrogate null distribution
permutate_behav <- function(data, cols_id) {
  data_ids <- unique(data[cols_id])
  data_ids_perm <- data_ids[sample.int(nrow(data_ids)), ]
  suff_tmp <- "_perm"
  names(data_ids_perm) <- paste0(cols_id, suff_tmp)
  bind_cols(data_ids, data_ids_perm) |>
    left_join(data, by = cols_id) |>
    select(-all_of(cols_id)) |>
    rename_with(
      ~ str_remove(.x, suff_tmp),
      ends_with(suff_tmp)
    )
}

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
        ~ mean(select(.x, {{ col_stats }}) > {{ col_stats }})
      ),
      .keep = "unused"
    )
}

extract_stats_cluster <- function(stats, .by,
                                  col_p_value = p.value,
                                  col_statistic = statistic,
                                  col_window = window_id,
                                  alpha = 0.05) {
  stats |>
    # order is essential for cluster detection
    arrange({{ col_window }}) |>
    summarise(
      find_largest_cluster({{ col_statistic }}, {{ col_p_value }} < alpha),
      .by = {{ .by }}
    )
}

find_largest_cluster <- function(statistic, is_sig) {
  clusters <- as_tibble(find_cluster(is_sig))
  if (nrow(clusters) == 0) {
    return(tibble(start = NA, end = NA, sum_t = 0))
  }
  clusters |>
    mutate(
      sum_t = map2_dbl(
        start, end,
        \(start, end) {
          sum(statistic[start:end])
        }
      )
    ) |>
    slice_max(sum_t)
}

find_cluster <- function(x, values_keep = 1) {
  # https://stackoverflow.com/a/43875717/5996475
  rle_x <- rle(x)
  end <- cumsum(rle_x$lengths)
  start <- c(1, lag(end)[-1] + 1)
  list(
    start = start[rle_x$values == values_keep],
    end = end[rle_x$values == values_keep]
  )
}
