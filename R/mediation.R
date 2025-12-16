fit_med <- function(model, data, ..., X = "X", Y = "Y", M = "M") {
  lavaan::sem(str_glue(str_c(model, collapse = "\n")), data)
}

combine_ccas <- function(data) {
  # pivot to wide for model fitting (creates igs_1, igs_2, ... and iss_1, ...)
  data_wider <- data |>
    pivot_wider(
      names_from = cca_id,
      values_from = c(igs, iss),
      names_glue = "{.value}_{cca_id}"
    )

  # helper: fit lm, extract standardized coefficients, and tidy them to (var, cca_id, coef)
  fit_weights <- function(fit) {
    parameters::model_parameters(fit, standardize = "refit") |>
      as.data.frame() |>
      filter(Parameter != "(Intercept)") |>
      separate_wider_delim(Parameter, "_", names = c("var", "cca_id")) |>
      mutate(cca_id = as.integer(cca_id)) |>
      select(var, cca_id, coef = Coefficient)
  }

  # get weights for iss and igs, then reshape to one row per cca_id
  weights <- bind_rows(
    fit_weights(lm(dprime ~ iss_1 + iss_2 + iss_3, data_wider)),
    fit_weights(lm(dprime ~ igs_1 + igs_2 + igs_3, data_wider))
  ) |>
    pivot_wider(
      names_from = var,
      values_from = coef,
      names_prefix = "weight_"
    )

  # join weights back and compute weighted summaries per subj_id / dprime
  data |>
    left_join(weights, by = "cca_id") |>
    summarise(
      c_igs = sum(igs * weight_igs) / sum(weight_igs),
      c_iss = sum(iss * weight_iss) / sum(weight_iss),
      .by = c(subj_id, dprime)
    )
}
