calc_iss <- function(patterns, pattern_semantics) {
  patterns |>
    mutate(
      iss = map_dbl(
        pattern,
        \(pattern) {
          cor(atanh(pattern), pattern_semantics, use = "pairwise")
        }
      ),
      .keep = "unused"
    )
}

calc_iss_stats <- function(data, ..., .by = c(cca_id, time_id)) {
  data |>
    summarise(
      broom::tidy(t.test(iss, ...)),
      .by = {{ .by }}
    )
}
