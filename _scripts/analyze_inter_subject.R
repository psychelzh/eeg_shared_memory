library(targets)

tar_option_set(
  packages = c("tidyverse"),
  controller = crew::crew_controller_local(
    name = "local",
    workers = 8
  ),
  garbage_collection = TRUE,
  memory = "transient",
  retrieval = "worker",
  storage = "worker"
)

calc_iss <- function(patterns_cca, pattern_semantics) {
  patterns_cca |>
    mutate(
      iss = map_dbl(
        pattern_cca,
        \(pattern_cca) {
          cor(atanh(pattern_cca), pattern_semantics, use = "pairwise")
        }
      ),
      .keep = "unused"
    )
}

calc_iss_stats <- function(data) {
  data |>
    summarise(
      broom::tidy(t.test(iss, alternative = "greater")),
      .by = c(cca_id, time_id)
    )
}

find_cluster <- function(statistic, p.value,
                         index = NULL,
                         keep = c("all", "largest"),
                         alpha = 0.05) {
  keep <- match.arg(keep)
  # https://stackoverflow.com/a/43875717/5996475
  rle_signif <- rle(p.value < alpha)
  if (!any(rle_signif$values)) {
    return(tibble(start = NA, end = NA, cluster_mass = 0))
  }
  end <- cumsum(rle_signif$lengths)
  start <- c(1, lag(end)[-1] + 1)
  clusters <- tibble(
    start = start[rle_signif$values],
    end = end[rle_signif$values]
  ) |>
    mutate(
      cluster_mass = map2_dbl(
        start, end,
        \(start, end) {
          sum(statistic[start:end])
        }
      )
    )
  if (!is.null(index)) {
    clusters$start <- index[clusters$start]
    clusters$end <- index[clusters$end]
  }
  if (keep == "largest") {
    clusters <- slice_max(clusters, cluster_mass)
  }
  clusters
}

list(
  tar_target(
    file_cca_y,
    "data/CorCAExtra/cca_y_subjs206_flat.parquet",
    format = "file"
  ),
  tar_target(
    subj_id_loop,
    arrow::open_dataset(file_cca_y) |>
      distinct(subj_id) |>
      pull(subj_id, as_vector = TRUE)
  ),
  tar_target(
    patterns_cca,
    arrow::open_dataset(file_cca_y) |>
      filter(subj_id == subj_id_loop) |>
      collect() |>
      pivot_wider(names_from = trial_id, values_from = y) |>
      reframe(
        pick(!time_id) |>
          slider::slide(
            \(x) as.dist(cor(x, use = "pairwise")),
            .before = 25,
            .after = 25,
            .step = 5,
            .complete = TRUE
          ) |>
          enframe(name = "time_id", value = "pattern_cca") |>
          filter(!map_lgl(pattern_cca, is.null)),
        .by = c(subj_id, cca_id)
      ),
    pattern = map(subj_id_loop)
  ),
  tar_target(file_seq, "config/sem_sequence.mat", format = "file"),
  tar_target(file_w2v, "data/stimuli/words_w2v.txt", format = "file"),
  tar_target(
    pattern_semantics,
    raveio::read_mat(file_seq)$SM[, 1:2] |>
      as_tibble(.name_repair = ~ c("trial_id", "word_id")) |>
      left_join(
        read_table(file_w2v, show_col_types = FALSE, col_names = FALSE),
        by = c("word_id" = "X1")
      ) |>
      filter(trial_id > 0) |>
      select(-word_id, -X2) |>
      column_to_rownames("trial_id") |>
      proxy::simil(method = "cosine")
  ),
  tar_target(iss, calc_iss(patterns_cca, pattern_semantics)),
  tar_target(iss_stats, calc_iss_stats(iss)),
  tarchetypes::tar_rep(
    iss_stats_permuted,
    calc_iss_stats(
      calc_iss(
        patterns_cca,
        seriation::permute(pattern_semantics, sample.int(150L))
      )
    ),
    reps = 10,
    batches = 100
  )
)
