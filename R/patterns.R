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

# intra/inter subject synchronizations ----
calc_sync_whole <- function(data) {
  data |>
    pivot_wider(names_from = subj_id, values_from = y_avg) |>
    summarise(
      pattern = list(calc_pattern(pick(matches("^\\d+$")))),
      .by = cca_id
    )
}

calc_sync_dynamic <- function(data) {
  data |>
    pivot_wider(names_from = subj_id, values_from = y_avg) |>
    reframe(
      calc_slide_window(pick(!time_id)),
      .by = cca_id
    )
}

calc_sync_between_halves <- function(data) {
  data |>
    pivot_wider(names_from = half, values_from = y_avg) |>
    reframe(
      {
        first <- pick(subj_id, time_id, first) |>
          pivot_wider(names_from = subj_id, values_from = first) |>
          column_to_rownames("time_id")
        second <- pick(subj_id, time_id, second) |>
          pivot_wider(names_from = subj_id, values_from = second) |>
          column_to_rownames("time_id")
        cor(first, second, use = "pairwise") |>
          atanh() |>
          as_tibble(rownames = "subj_id_first") |>
          pivot_longer(
            cols = -subj_id_first,
            names_to = "subj_id_second",
            values_to = "r"
          ) |>
          mutate(across(starts_with("subj_id"), as.integer))
      },
      .by = cca_id
    )
}

calc_sync_within_halves <- function(data) {
  data |>
    pivot_wider(names_from = subj_id, values_from = y_avg) |>
    reframe(
      pick(matches("^\\d+$")) |>
        cor(use = "pairwise") |>
        atanh() |>
        as_tibble(rownames = "subj_id_row") |>
        pivot_longer(cols = -subj_id_row, names_to = "subj_id_col", values_to = "r") |>
        mutate(across(starts_with("subj_id"), as.integer)) |>
        filter(subj_id_row < subj_id_col),
      .by = c(cca_id, half)
    )
}
