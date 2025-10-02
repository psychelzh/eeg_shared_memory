read_events_retrieval <- function(file, subjs) {
  read_tsv(file, show_col_types = FALSE) |>
    match_subj_id(subjs) |>
    # drop trials without response
    filter(resp != 0) |>
    # correct memory score
    mutate(score = as.double(xor(old_new == 1, resp >= 3)))
}

read_events_encoding <- function(file, subjs) {
  read_tsv(file, show_col_types = FALSE) |>
    match_subj_id(subjs) |>
    filter(word_id <= 150)
}

match_subj_id <- function(events, subjs) {
  events |>
    mutate(subj_id = match(subj, subjs)) |>
    filter(!is.na(subj_id))
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
    summarise(
      pc = mean(score == 1),
      .by = c(trial_id, word_id)
    )
}

calc_mem_content <- function(events, memorability) {
  events |>
    left_join(memorability, by = c("trial_id", "word_id")) |>
    # we should use old words only
    filter(word_id <= 150) |>
    mutate(resp = -resp) |>
    summarise(
      r = psych::polyserial(pick(pc), pick(resp))[, 1],
      .by = subj
    )
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
