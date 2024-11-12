calc_mediation <- function(igs, iss, perf) {
  igs |>
    inner_join(iss, by = c("subj_id", "cca_id")) |>
    inner_join(perf, by = "subj_id") |>
    nest(.by = cca_id) |>
    mutate(model = map(data, mediate_sem))
}
