calc_igs <- function(patterns_indiv, patterns_group) {
  patterns_indiv |>
    inner_join(patterns_group, by = c("subj_id", "cca_id")) |>
    mutate(
      igs = map2_dbl(pattern.x, pattern.y, cor, use = "pairwise"),
      .keep = "unused"
    )
}

calc_group_pattern <- function(data) {
  data |>
    summarise(
      y_avg = mean(y, na.rm = TRUE),
      .by = c(cca_id, trial_id, time_id)
    ) |>
    pivot_wider(
      names_from = trial_id,
      values_from = y_avg
    ) |>
    summarise(
      pattern = list(as.dist(cor(pick(matches("\\d+")), use = "pairwise"))),
      .by = cca_id
    )
}
