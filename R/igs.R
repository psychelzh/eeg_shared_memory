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
