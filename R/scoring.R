# memeory ability (performance) ----
calc_mem_perf <- function(events_retrieval) {
  count_trials <- events_retrieval |>
    distinct(word_id, old_new) |>
    count(old_new, name = "n_total")
  events_clean <- events_retrieval |>
    mutate(
      old_new = factor(old_new, c("old", "new")),
      response_type = factor(
        response_type,
        c("remember", "know", "unsure", "new")
      )
    ) |>
    filter(memory_type > 0)
  dprime_precise <- calc_dprime(
    events_clean, count_trials,
    types_signal = c("remember", "know")
  )
  dprime_adj <- dprime_precise |>
    summarise(
      across(
        c(hr, far),
        ~ .x[response_type == "know"] /
          (1 - .x[response_type == "remember"])
      ),
      .by = subj_id
    ) |>
    mutate(
      response_type = "knowadj",
      dprime = qnorm(hr) - qnorm(far)
    )
  dprime_avg <- bind_rows(dprime_precise, dprime_adj) |>
    summarise(
      dprime = mean(dprime[response_type != "know"]),
      .by = subj_id
    ) |>
    mutate(response_type = "avg_rk")
  dprime_coarse <- events_clean |>
    mutate(
      response_type = fct_collapse(
        response_type,
        old = c("remember", "know"),
        new = c("unsure", "new")
      )
    ) |>
    calc_dprime(count_trials, types_signal = "old")
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
  bind_rows(
    dprime_precise,
    dprime_adj,
    dprime_avg,
    dprime_coarse
  ) |>
    select(subj_id, index_name = response_type, score = dprime) |>
    bind_rows(grades)
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

# helper functions ----
calc_dprime <- function(data, count_trials, types_signal) {
  data |>
    count(subj_id, response_type, old_new, .drop = FALSE) |>
    left_join(count_trials, by = "old_new") |>
    mutate(
      rate = (n + 0.5) / (n_total + 1),
      .by = c(subj_id, response_type)
    ) |>
    filter(response_type %in% types_signal) |>
    mutate(
      type = case_match(
        old_new,
        "old" ~ "hr",
        "new" ~ "far"
      )
    ) |>
    pivot_wider(
      id_cols = c(subj_id, response_type),
      names_from = type,
      values_from = rate
    ) |>
    mutate(dprime = qnorm(hr) - qnorm(far))
}
