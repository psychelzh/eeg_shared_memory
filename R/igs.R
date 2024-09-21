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

calc_group_pattern <- function(data) {
  data |>
    summarise(
      y_avg = mean(y, na.rm = TRUE),
      .by = c(cca_id, trial_id, time_id)
    ) |>
    arrange(trial_id) |>
    pivot_wider(
      names_from = trial_id,
      values_from = y_avg
    ) |>
    summarise(
      pattern = list(
        atanh(as.dist(cor(pick(matches("\\d+")), use = "pairwise")))
      ),
      .by = cca_id
    )
}
