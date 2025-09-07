# group-averaged representations ----
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
      pattern = list(calc_pattern(pick(matches("\\d+")))),
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
      calc_slide_window(pick(!time_id)),
      .by = cca_id
    )
}

# individualized representations ----
calc_indiv_pattern <- function(data) {
  data |>
    pivot_wider(
      names_from = trial_id,
      values_from = y,
      names_sort = TRUE
    ) |>
    summarise(
      pattern = list(calc_pattern(pick(matches("^\\d+$")))),
      .by = c(subj_id, cca_id)
    )
}

calc_indiv_pattern_region <- function(data) {
  sz <- dim(data) # c(channels, time, trials, subjects)
  array(
    data,
    c(sz[1] * sz[2], sz[3], sz[4])
  ) |>
    apply(3, calc_pattern, simplify = FALSE) |>
    enframe(name = "subj_id", value = "pattern")
}

calc_indiv_pattern_dynamic <- function(data) {
  data |>
    pivot_wider(
      names_from = trial_id,
      values_from = y,
      names_sort = TRUE
    ) |>
    reframe(
      calc_slide_window(pick(!time_id)),
      .by = c(subj_id, cca_id)
    )
}

calc_indiv_pattern_dynamic_region <- function(data) {
  sz <- dim(data) # c(channels, time, trials, subjects)
  slider::slide(
    seq_len(sz[2]),
    \(idx) {
      array(
        data[, idx, , ],
        c(sz[1] * length(idx), sz[3], sz[4])
      ) |>
        apply(3, calc_pattern, simplify = FALSE) |>
        enframe(name = "subj_id", value = "pattern")
    },
    .before = 25,
    .after = 25,
    .step = 5,
    .complete = TRUE
  ) |>
    enframe(name = "time_id", value = "data") |>
    unnest(data)
}
