# IGS ----
calc_igs <- function(patterns_indiv, patterns_group) {
  by <- setdiff(
    intersect(names(patterns_indiv), names(patterns_group)),
    "pattern"
  )
  patterns_indiv |>
    inner_join(patterns_group, by = by) |>
    mutate(
      igs = atanh(map2_dbl(pattern.x, pattern.y, cor, use = "pairwise")),
      .keep = "unused"
    )
}

calc_igs_mem <- function(data_igs, mem_perf, ...) {
  correlate_mem_perf(data_igs, mem_perf, igs, ...)
}

fit_mem_pred <- function(mem_perf, ...) {
  data <- lst(...) |> # use `lst()` to keep argument names
    lapply(\(dat) rename(dat, x = last_col())) |> # data is in the last column
    bind_rows(.id = "src") |>
    pivot_wider(
      names_from = c(src, cca_id),
      values_from = x
    ) |>
    left_join(mem_perf, by = "subj_id") |>
    select(-subj_id)
  caret::train(
    dprime ~ .,
    data = data,
    method = "lm",
    trControl = caret::trainControl(method = "LOOCV")
  )
}

# ISS ----
calc_iss <- function(patterns_indiv, pattern_semantics) {
  patterns_indiv |>
    mutate(
      iss = pattern |>
        map_dbl(
          \(x) atanh(cor(x, pattern_semantics, use = "pairwise"))
        ) |>
        atanh(),
      .keep = "unused"
    )
}

# two data frames situation
calc_iss2 <- function(patterns_indiv, patterns_semantics) {
  by <- setdiff(
    intersect(names(patterns_indiv), names(patterns_semantics)),
    "pattern"
  )
  patterns_indiv |>
    inner_join(patterns_semantics, by = by) |>
    mutate(
      iss = atanh(map2_dbl(pattern.x, pattern.y, cor, use = "pairwise")),
      .keep = "unused"
    )
}

calc_iss_ind_r2 <- function(patterns_indiv, pattern_semantics) {
  patterns_indiv |>
    summarise(
      patterns = list(do.call(cbind, pattern)),
      .by = subj_id
    ) |>
    mutate(
      r2 = map_dbl(
        patterns,
        \(x) calc_r2_x_y(x, pattern_semantics)
      ),
      .keep = "unused"
    )
}

compare_iss <- function(data_iss) {
  data_iss |>
    mutate(cca_id = factor(cca_id)) |>
    lmerTest::lmer(iss ~ cca_id + (1 | subj_id), data = _) |>
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

calc_iss_mem <- function(data_iss, mem_perf, ...) {
  correlate_mem_perf(data_iss, mem_perf, iss, ...)
}

compare_iss_mem <- function(stats_iss_mem) {
  expand_grid(start = 1:3, end = 1:3) |>
    filter(start > end) |>
    mutate(
      map2(
        start,
        end,
        \(x, y) {
          with(
            stats_iss_mem,
            as_tibble(
              psych::r.test(206, estimate[[x]], estimate[[y]])[c("z", "p")]
            )
          )
        }
      ) |>
        list_rbind()
    )
}

# ISGS ----
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

# common functions ----
calc_stats_t <- function(data, col, ..., .by = c(cca_id, time_id)) {
  data |>
    summarise(
      broom::tidy(t.test({{ col }}, ...)),
      .by = {{ .by }}
    )
}
