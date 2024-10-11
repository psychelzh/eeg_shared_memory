calc_mem_perf <- function(data, subjs) {
  data |>
    mutate(acc = xor(old_new == 1, resp >= 3)) |>
    preproc.iquizoo:::calc_sdt(
      type_signal = 1,
      by = "subj",
      name_acc = "acc",
      name_type = "old_new"
    ) |>
    mutate(subj_id = match(subj, subjs)) |>
    filter(!is.na(subj_id)) |>
    select(subj_id, dprime)
}

calc_simil_mem <- function(mem_perf) {
  mem_perf |>
    column_to_rownames("subj_id") |>
    proxy::simil(method = "Euclidean")
}
