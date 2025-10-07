# statistics ----
calc_stats_t <- function(data, col, ..., .by = c(cca_id, time_id)) {
  data |>
    summarise(
      broom::tidy(t.test({{ col }}, ...)),
      .by = {{ .by }}
    )
}

regress_pattern <- function(y, x) {
  lm(as.vector(y) ~ as.vector(x), na.action = na.exclude) |>
    resid() |>
    vctrs::vec_restore(y)
}

calc_r2_x_y <- function(x, y) {
  rows_keep <- complete.cases(x, y)
  if (sum(rows_keep) < 2) {
    return(NA_real_)
  }
  x <- x[rows_keep, ]
  y <- y[rows_keep]
  cor(lm.fit(cbind(1, x), y)$fitted.values, y)^2
}

convert_p2_p1 <- function(
  p.value,
  statistic,
  alternative = c("greater", "less")
) {
  alternative <- match.arg(alternative)
  ifelse(
    xor(alternative == "greater", statistic > 0),
    1 - p.value / 2,
    p.value / 2
  )
}

compare_partial <- function(base, partial) {
  name_resp <- last(names(base))
  fit <- bind_rows(base = base, partial = partial, .id = "type") |>
    mutate(cca_id = factor(cca_id)) |>
    lmerTest::lmer(
      paste(name_resp, "~ type * cca_id + (1 | subj_id)"),
      data = _
    )
  emm <- emmeans::emmeans(fit, ~ type * cca_id)
  pairwise <- pairs(emm, simple = "type")
  list(
    stats = broom::tidy(emm),
    htest = broom::tidy(pairwise)
  )
}

# pattern related ----
calc_pattern <- function(x) {
  atanh(as.dist(cor(x, use = "pairwise")))
}

calc_inter_patterns <- function(pattern_x, pattern_y) {
  atanh(cor(pattern_x, pattern_y, use = "pairwise"))
}

calc_slide_window <- function(data, ...) {
  data |>
    slider::slide(
      calc_pattern,
      ...,
      .before = 25,
      .after = 25,
      .step = 5,
      .complete = TRUE
    ) |>
    enframe(name = "index", value = "pattern") |>
    filter(!map_lgl(pattern, is.null))
}

# misc ----
permute_dist <- function(dist) {
  seriation::permute(dist, sample.int(attr(dist, "Size")))
}

order_by_trial <- function(pattern, mapping) {
  order <- with(mapping, word_id[trial_id > 0])
  as.dist(pattern[order, order])
}

index_time_id <- function(time, onset = 51, sampling_rate = 256) {
  # here time is in seconds
  round(time * sampling_rate) + onset
}

calc_common_time_id_rt_locked <- function(
  medrts_encoding,
  time_id_range = c(1, 307)
) {
  range_time_id_rt <- range(medrts_encoding$time_id_rt)
  left_bound <- time_id_range[1] - range_time_id_rt[1]
  right_bound <- time_id_range[2] - range_time_id_rt[2]
  seq(left_bound, right_bound)
}

shift_time_id_rt_locked <- function(
  data,
  medrts_encoding,
  common_time_id_rt_locked
) {
  data |>
    inner_join(medrts_encoding, by = "subj_id") |>
    mutate(time_id = time_id - time_id_rt) |>
    filter(time_id %in% common_time_id_rt_locked) |>
    select(-time_id_rt, -medrt)
}
