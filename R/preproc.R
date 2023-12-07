calc_mem_perf <- function(events_retrieval) {
  count_trials <- events_retrieval |>
    distinct(word_id, old_new) |>
    count(old_new, name = "n_total")
  events_retrieval |>
    mutate(
      old_new = factor(old_new, c("old", "new")),
      response_type = factor(
        response_type,
        c("remember", "know", "unsure", "new")
      )
    ) |>
    filter(memory_type > 0) |>
    count(subj_id, response_type, old_new, .drop = FALSE) |>
    left_join(count_trials, by = "old_new") |>
    mutate(
      rate = (n + 0.5) / (n_total + 1),
      .by = c(subj_id, response_type)
    ) |>
    filter(response_type %in% c("remember", "know")) |>
    mutate(
      type = case_match(
        old_new,
        "old" ~ "hr",
        "new" ~ "far"
      )
    ) |>
    pivot_wider(
      id_cols = subj_id,
      names_from = c(response_type, type),
      values_from = rate
    ) |>
    mutate(
      dplyover::across2(
        contains("know"),
        contains("remember"),
        ~ .x / (1 - .y),
        .names = "knowadj_{suf}"
      )
    ) |>
    mutate(
      dplyover::across2(
        contains("hr"),
        contains("far"),
        ~ qnorm(.x) - qnorm(.y),
        .names = "dprime_{pre}"
      )
    ) |>
    select(subj_id, starts_with("dprime")) |>
    pivot_longer(
      starts_with("dprime"),
      names_to = c(".value", "mem_type"),
      names_pattern = "(.+)_(.+)"
    )
}

permutate_behav <- function(data, cols_id) {
  data_ids <- unique(data[cols_id])
  data_ids_perm <- data_ids[sample.int(nrow(data_ids)), ]
  suff_tmp <- "_perm"
  names(data_ids_perm) <- paste0(cols_id, suff_tmp)
  bind_cols(data_ids, data_ids_perm) |>
    left_join(data, by = cols_id) |>
    select(-all_of(cols_id)) |>
    rename_with(
      ~ str_remove(.x, suff_tmp),
      ends_with(suff_tmp)
    )
}
