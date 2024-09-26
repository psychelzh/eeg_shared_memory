calc_isps <- function(patterns_indiv) {
  patterns_indiv |>
    summarise(
      isps = list(calc_pattern(do.call(cbind, pattern))),
      .by = !c(subj_id, pattern)
    )
}
