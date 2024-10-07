calc_isps <- function(patterns_indiv) {
  patterns_indiv |>
    summarise(
      isps = list(calc_pattern(do.call(cbind, pattern))),
      .by = !c(subj_id, pattern)
    )
}

summarise_isps <- function(data_isps, se = FALSE) {
  data_isps |>
    mutate(
      isps_mean = map_dbl(isps, mean),
      isps_se = if (se) map_dbl(isps, \(x) sd(x) / sqrt(length(x))),
      .keep = "unused"
    )
}

calc_stats_isps <- function(summary_isps, summary_isps_permuted) {
  summary_isps_permuted |>
    rename(isps_mean_perm = isps_mean) |>
    select(!starts_with("tar")) |>
    chop(isps_mean_perm) |>
    left_join(summary_isps) |>
    mutate(
      p_perm = map2_dbl(
        isps_mean_perm, isps_mean,
        ~ (sum(.x >= .y) + 1) / (length(.x) + 1)
      )
    )
}
