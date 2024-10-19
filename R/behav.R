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

calc_mem_perf_precise <- function(data) {
  data |>
    count(subj_id, old_new, resp, .drop = FALSE) |>
    mutate(n = n + 0.5) |>
    mutate(
      rate = n / sum(n),
      .by = c(subj_id, old_new),
      .keep = "unused"
    ) |>
    filter(resp < 3) |>
    mutate(resp = factor(resp, labels = c("rem", "know"))) |>
    pivot_wider(names_from = resp, values_from = rate) |>
    # see Manns et al., 2003
    mutate(know = know / (1 - rem)) |>
    summarise(
      across(
        c(rem, know),
        \(x) qnorm(x[1]) - qnorm(x[2]),
        .names = "dprime_{.col}"
      ),
      .by = subj_id
    )
}
