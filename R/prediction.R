corr_mem <- function(
  data,
  mem_perf,
  ...,
  col_pred = last_col(),
  col_perf = dprime,
  by = "subj_id" # no need to support complex joins for now
) {
  name_pred <- names(select(data, {{ col_pred }}))
  name_perf <- names(select(mem_perf, {{ col_perf }}))
  data |>
    left_join(mem_perf, by = by) |>
    summarise(
      broom::tidy(cor.test(
        .data[[name_pred]],
        .data[[name_perf]],
        use = "pairwise",
        ...
      )),
      .by = !all_of(c(by, name_pred, name_perf))
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
