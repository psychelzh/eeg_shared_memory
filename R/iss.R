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
