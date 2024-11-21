fit_med <- function(model, data, ...,
                    X = "X", Y = "Y", M = "M") {
  lavaan::sem(str_glue(str_c(model, collapse = "\n")), data)
}

combine_data_ccas <- function(data) {
  data_weight <- data |>
    summarise(
      across(
        c(igs, iss),
        \(x) cor(x, dprime),
        .names = "weight_{.col}"
      ),
      .by = cca_id
    )
  data |>
    left_join(data_weight, by = "cca_id") |>
    summarise(
      igs = sum(igs * weight_igs) / sum(weight_igs),
      iss = sum(iss * weight_iss) / sum(weight_igs),
      .by = c(subj_id, dprime)
    )
}
