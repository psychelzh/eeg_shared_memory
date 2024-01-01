# memeory ability (performance) ----
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

calc_dist_mem_perf <- function(mem_perf) {
  mem_perf |>
    filter(mem_type %in% names(mem_types_report)) |>
    pivot_wider(
      names_from = mem_type,
      values_from = dprime
    ) |>
    column_to_rownames("subj_id") |>
    dist(method = "euclidean")
}

# similarity/distance of participants' responses ----
transform_resp_precise <- function(events_retrieval) {
  events_retrieval
}

transform_resp_coarse <- function(events_retrieval) {
  events_retrieval |>
    mutate(
      memory_type = case_match(
        memory_type,
        c(1, 2) ~ 1,
        c(3, 4) ~ 2,
        .default = 0
      )
    )
}

prepare_resp_mat <- function(resp, include) {
  if (include == "all") {
    include <- c("old", "new")
  }
  # note: 0's in memory_type are not removed now
  resp |>
    filter(old_new %in% include) |>
    pivot_wider(
      id_cols = subj_id,
      names_from = word_id,
      values_from = memory_type
    ) |>
    column_to_rownames("subj_id")
}

calc_dist_resp_mat <- function(resp_mat, method = c("sm", "gower")) {
  switch(method,
    sm = 1 - nomclust::sm(resp_mat),
    gower = resp_mat |>
      mutate(
        across(
          everything(),
          # 0 means no response and should be removed here
          \(x) factor(na_if(x, 0), ordered = TRUE)
        )
      ) |>
      proxy::simil(method = "Gower")
  )
}
