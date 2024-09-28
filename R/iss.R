calc_iss <- function(patterns_indiv, pattern_semantics) {
  patterns_indiv |>
    mutate(
      iss = pattern |>
        map_dbl(
          \(x) atanh(cor(x, pattern_semantics, use = "pairwise"))
        ) |>
        atanh(),
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

compare_iss <- function(data_iss) {
  data_iss |>
    mutate(cca_id = factor(cca_id)) |>
    lmerTest::lmer(iss ~ cca_id + (1 | subj_id), data = _) |>
    emmeans::emmeans(
      ~cca_id,
      lmer.df = "satterthwaite",
      lmerTest.limit = Inf
    ) |>
    emmeans::contrast("pairwise") |>
    broom::tidy() |>
    separate_wider_delim(
      contrast, " - ",
      names = c("start", "end")
    ) |>
    mutate(across(c("start", "end"), parse_number))
}

calc_iss_mem <- function(data_iss, mem_perf, ...) {
  correlate_mem_perf(data_iss, mem_perf, iss, ...)
}

compare_iss_mem <- function(stats_iss_mem) {
  expand_grid(start = 1:3, end = 1:3) |>
    filter(start > end) |>
    mutate(
      map2(
        start, end,
        \(x, y) {
          with(
            stats_iss_mem,
            as_tibble(
              psych::r.test(206, estimate[[x]], estimate[[y]])[c("z", "p")]
            )
          )
        }
      ) |>
        list_rbind()
    )
}
