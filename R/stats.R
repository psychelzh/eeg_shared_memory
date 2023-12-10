extract_stats_group <- function(file_parquet, mem_perf, ...) {
  arrow::open_dataset(file_parquet) |>
    left_join(mem_perf, by = "subj_id", relationship = "many-to-many") |>
    collect() |>
    summarise(
      cor.test(fisher_z, dprime, ...) |>
        broom::tidy() |>
        select(r = estimate, t = statistic, p = p.value),
      .by = c(region_id, trial_id, window_id, mem_type)
    )
}

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
