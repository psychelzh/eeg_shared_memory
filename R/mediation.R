calc_mediation <- function(model, data, ...) {
  data |>
    nest(.by = cca_id, .key = "data_cca") |>
    mutate(
      fit = map(
        data_cca,
        \(x) fit_med(model, x, ...)
      ),
      .keep = "unused"
    )
}

fit_med <- function(model, data, ...,
                    X = "X", Y = "Y", M = "M") {
  lavaan::sem(str_glue(str_c(model, collapse = "\n")), data)
}
