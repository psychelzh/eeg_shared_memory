convert_p2_p1 <- function(statistic, p.value,
                          alternative = c("greater", "less")) {
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
