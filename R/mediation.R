fit_med <- function(model, data, ..., X = "X", Y = "Y", M = "M") {
  lavaan::sem(str_glue(str_c(model, collapse = "\n")), data)
}

combine_data_ccas <- function(data) {
  data_wider <- data |>
    pivot_wider(names_from = cca_id, values_from = c(igs, iss))
  params_iss <- lm(dprime ~ iss_1 + iss_2 + iss_3, data_wider) |>
    parameters::model_parameters(standardize = "refit")
  params_igs <- lm(dprime ~ igs_1 + igs_2 + igs_3, data_wider) |>
    parameters::model_parameters(standardize = "refit")
  data_weight <- bind_rows(params_iss, params_igs) |>
    filter(Parameter != "(Intercept)") |>
    separate_wider_delim(Parameter, "_", names = c("var", "cca_id")) |>
    select(var, cca_id, coef = Coefficient) |>
    mutate(cca_id = as.integer(cca_id)) |>
    pivot_wider(names_from = var, values_from = coef, names_prefix = "weight_")
  data |>
    left_join(data_weight, by = "cca_id") |>
    summarise(
      igs = sum(igs * weight_igs) / sum(weight_igs),
      iss = sum(iss * weight_iss) / sum(weight_igs),
      .by = c(subj_id, dprime)
    )
}
