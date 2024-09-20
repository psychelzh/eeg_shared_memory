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

# for whole times series analysis, we would remove the first 200 ms baseline
index_onset <- floor(256 * (200 / 1000))

calc_iss <- function(patterns, pattern_semantics) {
  patterns |>
    mutate(
      iss = map_dbl(
        pattern,
        \(pattern) {
          cor(atanh(pattern), pattern_semantics, use = "pairwise")
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

calc_sync_smc <- function(sync_whole_trials, smc) {
  sync_whole_trials |>
    mutate(
      mantel = map(neu_sync, ~ vegan::mantel(.x, smc)),
      .keep = "unused"
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

tidy_mantel <- function(mantel) {
  tibble(
    statistic = mantel$statistic,
    p.value = mantel$signif,
    method = mantel$method
  )
}

list(
  tar_target(
    file_cca_y,
    "data/CorCAExtra/cca_y_subjs206.parquet",
    format = "file"
  ),
  tar_target(
    subj_id_loop,
    arrow::open_dataset(file_cca_y) |>
      distinct(subj_id) |>
      pull(subj_id, as_vector = TRUE)
  ),
  tar_target(
    patterns_indiv_dynamic,
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
          enframe(name = "time_id", value = "pattern") |>
          filter(!map_lgl(pattern, is.null)),
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
  tar_target(data_iss_dynamic, calc_iss(patterns_indiv_dynamic, pattern_semantics)),
  tar_target(stats_iss_dynamic, calc_iss_stats(data_iss_dynamic)),
  tarchetypes::tar_rep(
    data_iss_dynamic_permuted,
    calc_iss(
      patterns_indiv_dynamic,
      seriation::permute(pattern_semantics, sample.int(150L))
    ),
    reps = 10,
    batches = 100,
    iteration = "list"
  ),
  tarchetypes::tar_rep2(
    stats_iss_dynamic_permuted,
    calc_iss_stats(data_iss_dynamic_permuted, alternative = "greater"),
    data_iss_dynamic_permuted
  ),
  tar_target(
    clusters_stats_iss,
    stats_iss_dynamic |>
      mutate(p.value = convert_p2_p1(statistic, p.value)) |>
      calc_clusters_stats(stats_iss_dynamic_permuted)
  ),
  tar_target(
    patterns_indiv_whole,
    arrow::open_dataset(file_cca_y) |>
      filter(time_id >= index_onset) |>
      collect() |>
      pivot_wider(names_from = trial_id, values_from = y) |>
      summarise(
        pattern = list(as.dist(cor(pick(matches("^\\d+$")), use = "pairwise"))),
        .by = c(subj_id, cca_id)
      )
  ),
  tar_target(data_iss_whole, calc_iss(patterns_indiv_whole, pattern_semantics)),
  tar_target(stats_iss_whole, calc_iss_stats(data_iss_whole, .by = cca_id)),
  tar_target(
    iss_comparison,
    data_iss_whole |>
      mutate(cca_id = factor(cca_id)) |>
      lmerTest::lmer(iss ~ cca_id + (1 | subj_id), data = _) |>
      emmeans::emmeans(
        ~cca_id,
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
  ),
  tarchetypes::tar_file_read(
    subjs,
    "data/subj_206.txt",
    read = scan(!!.x)
  ),
  tarchetypes::tar_file_read(
    mem_perf,
    "data/behav/retrieval.tsv",
    read = read_tsv(!!.x, show_col_types = FALSE) |>
      mutate(acc = xor(old_new == 1, resp >= 3)) |>
      preproc.iquizoo:::calc_sdt(
        type_signal = 1,
        by = "subj",
        name_acc = "acc",
        name_type = "old_new"
      ) |>
      mutate(subj_id = match(subj, subjs)) |>
      filter(!is.na(subj_id)) |>
      select(subj_id, dprime)
  ),
  tar_target(
    stats_iss_mem_whole,
    data_iss_whole |>
      left_join(mem_perf, by = "subj_id") |>
      summarise(broom::tidy(cor.test(atanh(iss), dprime)), .by = cca_id)
  ),
  tar_target(
    comparison_iss_mem,
    expand_grid(start = 1:3, end = 1:3) |>
      filter(start > end) |>
      mutate(
        map2(
          start, end,
          \(x, y) {
            with(
              stats_iss_mem_whole,
              as_tibble(
                psych::r.test(206, estimate[[x]], estimate[[y]])[c("z", "p")]
              )
            )
          }
        ) |>
          list_rbind()
      )
  ),
  tar_target(
    stats_iss_mem_dynamic,
    data_iss_dynamic |>
      left_join(mem_perf, by = "subj_id") |>
      summarise(
        broom::tidy(cor.test(iss, dprime, use = "pairwise")),
        .by = c(cca_id, time_id)
      )
  ),
  tarchetypes::tar_rep(
    stats_iss_mem_dynamic_permuted,
    data_iss_dynamic |>
      left_join(
        mem_perf |>
          mutate(subj_id = sample(subj_id)),
        by = "subj_id"
      ) |>
      summarise(
        cor.test(iss, dprime, alternative = "greater", use = "pairwise") |>
          broom::tidy(),
        .by = c(cca_id, time_id)
      ),
    reps = 10,
    batches = 100
  ),
  tar_target(
    clusters_stats_iss_mem,
    stats_iss_mem_dynamic |>
      mutate(p.value = convert_p2_p1(statistic, p.value)) |>
      calc_clusters_stats(stats_iss_mem_dynamic_permuted)
  ),
  tar_target(
    cca_y_halves_trials,
    arrow::open_dataset(file_cca_y) |>
      mutate(half = if_else(trial_id <= 75, "first", "second")) |>
      filter(!is.nan(y)) |>
      count(subj_id, cca_id, time_id, half) |>
      distinct(subj_id, cca_id, half, n) |>
      collect()
  ),
  tar_target(
    cca_y_halves,
    arrow::open_dataset(file_cca_y) |>
      mutate(half = if_else(trial_id <= 75, "first", "second")) |>
      filter(!is.nan(y)) |>
      summarise(
        y_avg = mean(y),
        .by = c(subj_id, cca_id, time_id, half)
      ) |>
      collect()
  ),
  tar_target(
    sync_inter_subjs,
    cca_y_halves |>
      filter(time_id >= index_onset) |>
      pivot_wider(names_from = subj_id, values_from = y_avg) |>
      reframe(
        cor(pick(matches("^\\d+$")), use = "pairwise") |>
          as_tibble(rownames = "row") |>
          pivot_longer(cols = -row, names_to = "col", values_to = "r") |>
          mutate(across(c(row, col), as.integer)) |>
          filter(row < col),
        .by = c(cca_id, half)
      )
  ),
  tar_target(
    sync_inter_halves,
    cca_y_halves |>
      filter(time_id >= index_onset) |>
      pivot_wider(names_from = half, values_from = y_avg) |>
      reframe(
        {
          first <- pick(subj_id, time_id, first) |>
            pivot_wider(names_from = subj_id, values_from = first) |>
            column_to_rownames("time_id")
          second <- pick(subj_id, time_id, second) |>
            pivot_wider(names_from = subj_id, values_from = second) |>
            column_to_rownames("time_id")
          cor(first, second, use = "pairwise") |>
            as_tibble(rownames = "first") |>
            pivot_longer(cols = -first, names_to = "second", values_to = "r") |>
            mutate(across(c(first, second), as.integer))
        },
        .by = cca_id
      )
  ),
  tar_target(
    whole_erps,
    arrow::read_parquet(file_cca_y) |>
      summarise(
        y_avg = mean(y, na.rm = TRUE),
        .by = c(subj_id, cca_id, time_id)
      )
  ),
  tarchetypes::tar_file_read(
    smc,
    "data/behav/simil.rds",
    read = readRDS(!!.x)$mat[[4]]
  ),
  tar_target(
    sync_whole_trials,
    whole_erps |>
      filter(time_id >= index_onset) |>
      pivot_wider(names_from = subj_id, values_from = y_avg) |>
      summarise(
        neu_sync = list(cor(pick(matches("^\\d+$")), use = "pairwise")),
        .by = cca_id
      )
  ),
  tar_target(
    sync_smc,
    calc_sync_smc(sync_whole_trials, smc)
  ),
  tar_target(
    sync_dynamic,
    whole_erps |>
      pivot_wider(names_from = subj_id, values_from = y_avg) |>
      reframe(
        pick(!time_id) |>
          slider::slide(
            \(x) as.dist(cor(x, use = "pairwise")),
            .before = 25,
            .after = 25,
            .step = 5,
            .complete = TRUE
          ) |>
          enframe(name = "time_id", value = "neu_sync") |>
          filter(!map_lgl(neu_sync, is.null)),
        .by = cca_id
      )
  ),
  tar_target(
    sync_smc_dynamic,
    calc_sync_smc(sync_dynamic, smc)
  ),
  tarchetypes::tar_rep(
    sync_smc_dynamic_permuted,
    calc_sync_smc(
      sync_dynamic,
      seriation::permute(smc, sample.int(206L))
    ),
    reps = 10,
    batches = 100
  ),
  tarchetypes::tar_rep2(
    stats_sync_smc_dynamic_permuted,
    sync_smc_dynamic_permuted |>
      mutate(
        map(mantel, tidy_mantel) |>
          list_rbind(),
        .keep = "unused"
      ),
    sync_smc_dynamic_permuted
  ),
  tar_target(
    stats_sync_smc_dynamic,
    sync_smc_dynamic |>
      mutate(
        map(mantel, tidy_mantel) |>
          list_rbind(),
        .keep = "unused"
      )
  ),
  tar_target(
    clusters_stats_sync_smc_dynamic,
    stats_sync_smc_dynamic |>
      calc_clusters_stats(stats_sync_smc_dynamic_permuted)
  )
)
