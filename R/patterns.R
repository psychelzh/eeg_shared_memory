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

regress_patterns <- function(patterns_y, patterns_x, by = NULL) {
  if (is.null(by)) {
    by <- intersect(names(patterns_y), names(patterns_x)) |>
      setdiff("pattern")
  }
  patterns_x |>
    left_join(patterns_y, by = by) |>
    mutate(
      pattern = map2(pattern.y, pattern.x, regress_pattern),
      .keep = "unused"
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

prepare_sync_inter_intra <- function(sync_between_halves) {
  sync_between_halves |>
    mutate(
      type = case_when(
        subj_id_first == subj_id_second ~ "intra",
        subj_id_first < subj_id_second ~ "inter_ahead",
        TRUE ~ "inter_behind"
      ),
      subj_id = if_else(
        type != "inter_behind",
        subj_id_first,
        subj_id_second
      ),
      .keep = "unused"
    ) |>
    summarise(
      sync = mean(r),
      .by = c(cca_id, subj_id, type)
    ) |>
    # here we use "inter_ahead" only
    filter(type %in% c("intra", "inter_ahead")) |>
    mutate(
      cca_id = factor(cca_id),
      type = factor(type, levels = c("intra", "inter_ahead"))
    )
}
