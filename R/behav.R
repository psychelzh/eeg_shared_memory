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
    # Count directly so unobserved responses contribute zero, not missing rows.
    # Preserve four-category smoothing: 0.5 per response, hence N + 2.
    summarise(
      rem = (sum(resp == 1) + 0.5) / (n() + 2),
      know = (sum(resp == 2) + 0.5) / (n() + 2),
      .by = c(subj_id, old_new)
    ) |>
    # Independent Remember/Know adjustment; see Manns et al., 2003.
    mutate(know = know / (1 - rem)) |>
    pivot_wider(
      names_from = old_new,
      values_from = c(rem, know)
    ) |>
    # Stimulus codes: 1 = old (signal), 2 = new (noise).
    transmute(
      subj_id,
      dprime_rem = qnorm(rem_1) - qnorm(rem_2),
      dprime_know = qnorm(know_1) - qnorm(know_2)
    )
}

calc_medrt_encoding <- function(data) {
  data |>
    summarise(
      medrt = median(rt[acc == 1], na.rm = TRUE),
      .by = subj_id
    ) |>
    mutate(time_id_rt = index_time_id(medrt))
}
