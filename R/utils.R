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

get_resid <- function(y, x) {
  resid(lm(y ~ x, na.action = na.exclude))
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
    enframe(name = "time_id", value = "pattern") |>
    filter(!map_lgl(pattern, is.null))
}

calc_pattern <- function(x) {
  atanh(as.dist(cor(x, use = "pairwise")))
}

calc_inter_patterns <- function(pattern_x, pattern_y) {
  atanh(cor(pattern_x, pattern_y, use = "pairwise"))
}

permute_dist <- function(dist) {
  seriation::permute(dist, sample.int(attr(dist, "Size")))
}

regress_pattern <- function(y, x) {
  lm(as.vector(y) ~ as.vector(x), na.action = na.exclude) |>
    resid() |>
    vctrs::vec_restore(y)
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

calc_r2_x_y <- function(x, y) {
  rows_keep <- complete.cases(x, y)
  if (sum(rows_keep) < 2) {
    return(NA_real_)
  }
  x <- x[rows_keep, ]
  y <- y[rows_keep]
  cor(lm.fit(cbind(1, x), y)$fitted.values, y)^2
}

order_by_trial <- function(pattern, mapping) {
  order <- with(mapping, word_id[trial_id > 0])
  as.dist(pattern[order, order])
}
