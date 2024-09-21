convert_p2_p1 <- function(statistic, p.value,
                          alternative = c("greater", "less")) {
  alternative <- match.arg(alternative)
  ifelse(
    xor(alternative == "greater", statistic > 0),
    1 - p.value / 2,
    p.value / 2
  )
}

tidy_mantel <- function(mantel) {
  tibble(
    statistic = mantel$statistic,
    p.value = mantel$signif,
    method = mantel$method
  )
}

get_resid <- function(y, x) {
  resid(lm(y ~ x, na.action = na.exclude))
}
