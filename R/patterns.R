calc_group_pattern <- function(data) {
  data |>
    summarise(
      y_avg = mean(y, na.rm = TRUE),
      .by = c(cca_id, trial_id, time_id)
    ) |>
    pivot_wider(
      names_from = trial_id,
      values_from = y_avg,
      names_sort = TRUE
    ) |>
    summarise(
      pattern = list(.calc_pattern(pick(matches("\\d+")))),
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
      values_from = y_avg,
      names_sort = TRUE
    ) |>
    reframe(
      .calc_pattern_dynamic(pick(!time_id)),
      .by = cca_id
    )
}

calc_indiv_pattern <- function(data) {
  data |>
    pivot_wider(
      names_from = trial_id,
      values_from = y,
      names_sort = TRUE
    ) |>
    summarise(
      pattern = list(.calc_pattern(pick(matches("^\\d+$")))),
      .by = c(subj_id, cca_id)
    )
}

calc_indiv_pattern_dynamic <- function(data) {
  data |>
    pivot_wider(
      names_from = trial_id,
      values_from = y,
      names_sort = TRUE
    ) |>
    reframe(
      .calc_pattern_dynamic(pick(!time_id)),
      .by = c(subj_id, cca_id)
    )
}

# helper functions
.calc_pattern_dynamic <- function(data) {
  data |>
    slider::slide(
      .calc_pattern,
      .before = 25,
      .after = 25,
      .step = 5,
      .complete = TRUE
    ) |>
    enframe(name = "time_id", value = "pattern") |>
    filter(!map_lgl(pattern, is.null))
}

.calc_pattern <- function(x) {
  atanh(as.dist(cor(x, use = "pairwise")))
}
