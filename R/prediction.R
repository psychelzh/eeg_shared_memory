corr_mem <- function(
  data,
  mem_perf,
  ...,
  col_pred = last_col(),
  col_perf = last_col(),
  by = "subj_id" # no need to support complex joins for now
) {
  name_pred <- names(select(data, {{ col_pred }}))
  name_perf <- names(select(mem_perf, {{ col_perf }}))
  # assume `mem_perf` does not contain grouping variables
  # if so, split them out
  vars_grouping <- setdiff(names(data), c(by, name_pred))
  data |>
    left_join(mem_perf, by = by) |>
    summarise(
      broom::tidy(cor.test(
        .data[[name_pred]],
        .data[[name_perf]],
        use = "pairwise",
        ...
      )),
      .by = all_of(vars_grouping)
    )
}

compare_corr_mem <- function(stats, ..., n = 206) {
  expand_grid(start = 1:3, end = 1:3) |>
    filter(start > end) |>
    mutate(
      map2(
        start,
        end,
        \(x, y) {
          with(
            stats,
            psych::r.test(n, estimate[[x]], estimate[[y]])[c("z", "p")] |>
              as_tibble()
          )
        }
      ) |>
        list_rbind()
    )
}

fit_mem_pred <- function(mem_perf, ...) {
  data <- lst(...) |> # use `lst()` to keep argument names
    lapply(\(dat) rename(dat, x = last_col())) |> # data is in the last column
    bind_rows(.id = "src") |>
    pivot_wider(
      names_from = c(src, cca_id),
      values_from = x
    ) |>
    left_join(mem_perf, by = "subj_id") |>
    select(-subj_id)
  caret::train(
    dprime ~ .,
    data = data,
    method = "lm",
    trControl = caret::trainControl(method = "LOOCV")
  )
}
