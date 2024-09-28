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

compare_igs <- function(igs, igs_partial) {
  fit <- bind_rows(igs = igs, partial = igs_partial, .id = "type") |>
    mutate(cca_id = factor(cca_id)) |>
    lmerTest::lmer(igs ~ type * cca_id + (1 | subj_id), data = _)
  emm <- emmeans::emmeans(fit, ~ type * cca_id)
  pairwise <- pairs(emm, simple = "type")
  list(
    stats = broom::tidy(emm),
    htest = broom::tidy(pairwise)
  )
}

calc_igs_mem <- function(data_igs, mem_perf, ...) {
  correlate_mem_perf(data_igs, mem_perf, igs, ...)
}
