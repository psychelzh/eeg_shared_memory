library(targets)

tar_option_set(
  packages = c("tidyverse"),
  controller = crew::crew_controller_local(
    name = "local",
    workers = 12
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

calc_iss_stats <- function(data, ..., .by = c(cca_id, time_id)) {
  data |>
    summarise(
      broom::tidy(t.test(iss, ...)),
      .by = {{ .by }}
    )
}

calc_clusters_stats <- function(stats, stats_permuted,
                                by = "cca_id",
                                col_statistic = statistic,
                                col_p_value = p.value,
                                col_time_id = time_id,
                                col_id_permuted = starts_with("tar")) {
  clusters <- stats |>
    reframe(
      find_cluster(
        {{ col_statistic }},
        {{ col_p_value }},
        {{ col_time_id }}
      ),
      .by = all_of(by)
    )
  clusters_permuted <- stats_permuted |>
    reframe(
      find_cluster(
        {{ col_statistic }},
        {{ col_p_value }},
        {{ col_time_id }},
        keep = "largest"
      ),
      .by = c(all_of(by), {{ col_id_permuted }})
    )
  clusters |>
    left_join(
      clusters_permuted |>
        select(cca_id, cluster_mass_perm = cluster_mass) |>
        chop(cluster_mass_perm),
      by = by
    ) |>
    mutate(
      p_perm = map2_dbl(
        cluster_mass_perm,
        cluster_mass,
        ~ (sum(.x >= .y) + 1) / (length(.x) + 1)
      )
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

convert_p2_p1 <- function(statistic, p.value,
                          alternative = c("greater", "less")) {
  alternative <- match.arg(alternative)
  ifelse(
    xor(alternative == "greater", statistic > 0),
    1 - p.value / 2,
    p.value / 2
  )
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
  tar_target(data_iss, calc_iss(patterns_cca, pattern_semantics)),
  tar_target(iss_stats, calc_iss_stats(data_iss)),
  tarchetypes::tar_rep(
    data_iss_permuted,
    calc_iss(
      patterns_cca,
      seriation::permute(pattern_semantics, sample.int(150L))
    ),
    reps = 10,
    batches = 100,
    iteration = "list"
  ),
  tarchetypes::tar_rep2(
    iss_stats_permuted,
    calc_iss_stats(data_iss_permuted, alternative = "greater"),
    data_iss_permuted
  ),
  tar_target(
    clusters_stats_iss,
    iss_stats |>
      mutate(p.value = convert_p2_p1(statistic, p.value)) |>
      calc_clusters_stats(iss_stats_permuted)
  ),
  tar_target(
    patterns_cca_whole,
    arrow::open_dataset(file_cca_y) |>
      filter(time_id >= 51) |>
      collect() |>
      pivot_wider(names_from = trial_id, values_from = y) |>
      summarise(
        pattern_cca = list(as.dist(cor(pick(matches("^\\d+$"))))),
        .by = c(subj_id, cca_id)
      )
  ),
  tar_target(data_iss_whole, calc_iss(patterns_cca_whole, pattern_semantics)),
  tar_target(iss_stats_whole, calc_iss_stats(data_iss_whole, .by = cca_id)),
  tar_target(
    iss_comparison,
    data_iss_whole |>
      mutate(cca_id = factor(cca_id)) |>
      lmerTest::lmer(iss ~ cca_id + (1 | subj_id), data = _) |>
      emmeans::emmeans(
        ~ cca_id,
        lmer.df = "satterthwaite",
        lmerTest.limit = Inf
      ) |>
      emmeans::contrast("pairwise") |>
      broom::tidy() |>
      separate_wider_delim(
        contrast, " - ",
        names = c("start", "end")
      ) |>
      mutate(across(c("start", "end"), parse_number))
  )
)
