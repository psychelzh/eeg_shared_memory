clean_events <- function(file, subjs) {
  read_tsv(file, show_col_types = FALSE) |>
    mutate(subj_id = match(subj, subjs)) |>
    # drop trials without response
    filter(resp != 0, !is.na(subj_id)) |>
    # correct memory score
    mutate(score = as.double(xor(old_new == 1, resp >= 3)))
}

calc_mem_perf <- function(data) {
  data |>
    preproc.iquizoo:::calc_sdt(
      type_signal = 1,
      by = "subj_id",
      name_acc = "score",
      name_type = "old_new"
    ) |>
    select(subj_id, dprime)
}

calc_simil_mem <- function(mem_perf) {
  mem_perf |>
    column_to_rownames("subj_id") |>
    proxy::simil(method = "Euclidean")
}

calc_memorability <- function(data) {
  data |>
    summarise(pc = mean(score == 1), .by = trial_id)
}
