calc_igs <- function(patterns_indiv, patterns_group) {
  by <- setdiff(
    intersect(names(patterns_indiv), names(patterns_group)),
    "pattern"
  )
  patterns_indiv |>
    inner_join(patterns_group, by = by) |>
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

calc_group_pattern_dynamic <- function(data) {
  data |>
    summarise(
      y_avg = mean(y, na.rm = TRUE),
      .by = c(cca_id, trial_id, time_id)
    ) |>
    pivot_wider(
      names_from = trial_id,
      values_from = y_avg
    ) |>
    reframe(
      pick(!time_id) |>
        slider::slide(
          \(x) as.dist(cor(x, use = "pairwise")),
          .before = 25,
          .after = 25,
          .step = 5,
          .complete = TRUE
        ) |>
        enframe(name = "time_id", value = "pattern") |>
        filter(!map_lgl(pattern, is.null)),
      .by = cca_id
    )
}
