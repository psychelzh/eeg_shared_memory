# memeory ability (performance) ----
calc_mem_perf <- function(events_retrieval) {
  count_trials <- events_retrieval |>
    distinct(word_id, old_new) |>
    count(old_new, name = "n_total")
  dprimes  <- events_retrieval |>
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
        .names = "{pre}"
      )
    ) |>
    select(subj_id, !contains("_")) |>
    mutate(avg_rk = (remember + knowadj) / 2) |>
    pivot_longer(
      !subj_id,
      names_to = "index_name",
      values_to = "score"
    )
  grades <- events_retrieval |>
    filter(memory_type != 0) |>
    mutate(
      score = if_else(
        old_new == "old",
        5 - memory_type,
        memory_type
      )
    ) |>
    summarise(
      score = mean(score),
      .by = subj_id
    ) |>
    add_column(index_name = "avg_score", .before = "score")
  bind_rows(dprimes, grades)
}

calc_dist_mem_perf <- function(mem_perf, basis = c("knowadj", "remember")) {
  mat <- mem_perf |>
    filter(index_name %in% basis) |>
    pivot_wider(
      names_from = index_name,
      values_from = score
    ) |>
    column_to_rownames("subj_id")
  # use separate expression to make a clean "call" attributes
  dist(mat, method = "euclidean")
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
