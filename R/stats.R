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
