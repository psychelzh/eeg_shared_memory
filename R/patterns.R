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
  data_avg_wider <- data |>
    summarise(
      y_avg = mean(y, na.rm = TRUE),
      .by = c(cca_id, trial_id, time_id)
    ) |>
    pivot_wider(
      names_from = trial_id,
      values_from = y_avg,
      names_sort = TRUE
    )
  res <- data_avg_wider |>
    reframe(
      calc_slide_window(pick(!time_id)),
      .by = cca_id
    )
  data_avg_wider |>
    select(cca_id, time_id) |>
    # ensure the original time_id is kept
    mutate(
      index = row_number(time_id),
      .before = time_id,
      .by = cca_id
    ) |>
    inner_join(res, by = c("cca_id", "index")) |>
    select(-index)
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
  data_wider <- data |>
    pivot_wider(
      names_from = trial_id,
      values_from = y,
      names_sort = TRUE
    )
  res <- data_wider |>
    reframe(
      calc_slide_window(pick(!time_id)),
      .by = c(subj_id, cca_id)
    )
  data_wider |>
    select(subj_id, cca_id, time_id) |>
    # ensure the original time_id is kept
    mutate(
      index = row_number(time_id),
      .before = time_id,
      .by = c(subj_id, cca_id)
    ) |>
    inner_join(res, by = c("subj_id", "cca_id", "index")) |>
    select(-index)
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

# correlation between patterns ----
## single y pattern situation ----
corr_patterns_1 <- function(patterns_x, pattern_y, name) {
  patterns_x |>
    mutate(
      "{name}" := map_dbl(
        pattern,
        calc_inter_patterns,
        pattern_y
      ),
      .keep = "unused"
    )
}

## two data frames situation ----
corr_patterns <- function(patterns_x, patterns_y, name) {
  by <- setdiff(
    intersect(names(patterns_x), names(patterns_y)),
    "pattern"
  )
  if (length(by) == 0) {
    data <- cross_join(patterns_x, patterns_y)
  } else {
    data <- inner_join(patterns_x, patterns_y, by = by)
  }
  data |>
    mutate(
      "{name}" := map2_dbl(
        pattern.x,
        pattern.y,
        calc_inter_patterns
      ),
      .keep = "unused"
    )
}

## bind patterns among CCA components to calculate R2 ----
corr_patterns_r2 <- function(patterns_x, pattern_y) {
  patterns_x |>
    summarise(
      patterns = list(do.call(cbind, pattern)),
      .by = subj_id
    ) |>
    mutate(
      r2 = map_dbl(
        patterns,
        \(x) calc_r2_x_y(x, pattern_y)
      ),
      .keep = "unused"
    )
}

## comparisons ----
compare_corr_patterns <- function(data, col = last_col()) {
  name_resp <- names(select(data, {{ col }}))
  data |>
    mutate(cca_id = factor(cca_id)) |>
    lmerTest::lmer(
      str_c(name_resp, " ~ cca_id + (1 | subj_id)"),
      data = _
    ) |>
    emmeans::emmeans(
      ~cca_id,
      lmer.df = "satterthwaite",
      lmerTest.limit = Inf
    ) |>
    emmeans::contrast("pairwise") |>
    broom::tidy() |>
    separate_wider_delim(
      contrast,
      " - ",
      names = c("start", "end")
    ) |>
    mutate(across(c("start", "end"), parse_number))
}

# inter-subject pattern similarity ----
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
        isps_mean_perm,
        isps_mean,
        ~ (sum(.x >= .y) + 1) / (length(.x) + 1)
      )
    )
}

# regressions for control analysis ----
regress_patterns_1 <- function(patterns_y, pattern_x) {
  patterns_y |>
    mutate(pattern = map(pattern, regress_pattern, pattern_x))
}

regress_patterns_1r <- function(pattern_y, patterns_x) {
  patterns_x |>
    mutate(pattern = map(pattern, \(x) regress_pattern(pattern_y, x)))
}

regress_patterns <- function(patterns_y, patterns_x) {
  by <- setdiff(
    intersect(names(patterns_x), names(patterns_y)),
    "pattern"
  )
  if (length(by) == 0) {
    data <- cross_join(patterns_y, patterns_x)
  } else {
    data <- inner_join(patterns_y, patterns_x, by = by)
  }
  data |>
    mutate(
      pattern = map2(pattern.y, pattern.x, regress_pattern),
      .keep = "unused"
    )
}
