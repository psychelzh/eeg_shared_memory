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

calc_slide_window <- function(data, fun, name_value, ...) {
  data |>
    slider::slide(
      fun,
      ...,
      .before = 25,
      .after = 25,
      .step = 5,
      .complete = TRUE
    ) |>
    enframe(name = "time_id", value = name_value) |>
    filter(!map_lgl(.data[[name_value]], is.null))
}

calc_pattern <- function(x, fisher_z = TRUE) {
  ret <- as.dist(cor(x, use = "pairwise"))
  if (fisher_z) {
    atanh(ret)
  } else {
    ret
  }
}
