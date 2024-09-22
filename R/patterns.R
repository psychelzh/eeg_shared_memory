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
      .calc_pattern_dynamic(pick(!time_id)),
      .by = cca_id
    )
}

calc_indiv_pattern <- function(data) {
  data |>
    pivot_wider(names_from = trial_id, values_from = y) |>
    summarise(
      pattern = list(as.dist(cor(pick(matches("^\\d+$")), use = "pairwise"))),
      .by = c(subj_id, cca_id)
    )
}

calc_indiv_pattern_dynamic <- function(data) {
  data |>
    pivot_wider(names_from = trial_id, values_from = y) |>
    reframe(
      .calc_pattern_dynamic(pick(!time_id)),
      .by = c(subj_id, cca_id)
    )
}

.calc_pattern_dynamic <- function(data) {
  data |>
    slider::slide(
      \(x) as.dist(cor(x, use = "pairwise")),
      .before = 25,
      .after = 25,
      .step = 5,
      .complete = TRUE
    ) |>
    enframe(name = "time_id", value = "pattern") |>
    filter(!map_lgl(pattern, is.null))
}
