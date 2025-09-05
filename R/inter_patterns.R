# single y pattern situation ----
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

# two data frames situation ----
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

# ISPS ----
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

compare_inter_patterns <- function(data, col = last_col()) {
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
