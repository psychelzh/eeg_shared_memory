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


extract_response_shared <- function(events_encoding, events_retrieval) {
  events_retrieval |>
    filter(memory_type != 0) |>
    pivot_wider(
      id_cols = word_id,
      names_from = subj_id,
      values_from = response_type
    ) |>
    mutate(
      dplyover::across2x(
        -1, -1,
        ~ if_else(.x == .y, .x, NA),
        .comb = "minimal"
      )
    ) |>
    select(contains("_")) |>
    pivot_longer(
      -word_id,
      names_to = "subj_pair",
      values_to = "response_type_shared"
    ) |>
    separate(
      subj_pair,
      c("subj_id_col", "subj_id_row"),
      convert = TRUE
    ) |>
    mutate(
      response_type_shared = case_match(
        response_type_shared,
        "remember" ~ "Rem",
        "know" ~ "Know",
        "unsure" ~ "Unsure",
        "new" ~ "New",
        .ptype = factor(levels = c("Rem", "Know", "Unsure", "New"))
      )
    ) |>
    nest(.by = word_id, .key = "resp_matched") |>
    left_join(
      distinct(events_encoding, trial_id, word_id, word_category),
      by = "word_id"
    )
}

filter_shared <- function(file, response_shared) {
  arrow::read_parquet(file) |>
    inner_join(response_shared, by = "trial_id") |>
    mutate(
      filtered = map2(
        resp_matched,
        fisher_z,
        ~ .x |>
          add_column(fisher_z = .y) |>
          filter(!is.na(response_type_shared))
      ),
      .keep = "unused"
    ) |>
    unnest(filtered) |>
    filter(
      sum(!is.na(fisher_z)) >= 5,
      .by = c(region_id, word_category, contains("subj_id"))
    )
}
