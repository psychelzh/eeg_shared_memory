correlate_mem_perf <- function(data, mem_perf, col, ...) {
  data |>
    left_join(mem_perf, by = "subj_id") |>
    summarise(
      broom::tidy(cor.test({{ col }}, dprime, use = "pairwise", ...)),
      .by = !c(subj_id, {{ col }}, dprime)
    )
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

permute_dist <- function(dist) {
  seriation::permute(dist, sample.int(attr(dist, "Size")))
}

regress_pattern <- function(y, x) {
  lm(as.vector(y) ~ as.vector(x), na.action = na.exclude) |>
    resid() |>
    vctrs::vec_restore(y)
}

compare_partial <- function(base, partial) {
  fit <- bind_rows(base = base, partial = partial, .id = "type") |>
    mutate(cca_id = factor(cca_id)) |>
    lmerTest::lmer(igs ~ type * cca_id + (1 | subj_id), data = _)
  emm <- emmeans::emmeans(fit, ~ type * cca_id)
  pairwise <- pairs(emm, simple = "type")
  list(
    stats = broom::tidy(emm),
    htest = broom::tidy(pairwise)
  )
}
