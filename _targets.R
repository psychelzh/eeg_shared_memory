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

tar_source()

# for whole times series analysis, we would remove the first 200 ms baseline
index_onset <- floor(256 * (200 / 1000))
num_subj <- 206L
config_num_subjs <- data.frame(size = seq(20, num_subj, by = 20)) |>
  dplyr::mutate(paired = size <= 100)

targets_patterns_group_whole_resampled <- tarchetypes::tar_map(
  config_num_subjs,
  tarchetypes::tar_rep(
    subjs_sampled,
    resample(num_subj, size, paired),
    reps = 10,
    batches = 10,
    iteration = "list"
  ),
  tarchetypes::tar_rep2(
    patterns_group_whole_resampled,
    lapply(
      zutils::select_list(subjs_sampled, !starts_with("tar")),
      \(subjs) {
        arrow::read_parquet(file_cca_y) |>
          filter(time_id >= index_onset, subj_id %in% subjs) |>
          calc_group_pattern()
      }
    ) |>
      list_rbind(names_to = "pair"),
    subjs_sampled
  ),
  tarchetypes::tar_rep2(
    patterns_group_stability,
    if (paired) {
      patterns_group_whole_resampled |>
        pivot_wider(names_from = pair, values_from = pattern) |>
        mutate(r = map2_dbl(`1`, `2`, cor), .keep = "unused")
    } else {
      tibble()
    },
    patterns_group_whole_resampled
  )
)

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
  targets_patterns_group_whole_resampled,
  tarchetypes::tar_combine(
    patterns_group_stability,
    targets_patterns_group_whole_resampled$patterns_group_stability,
    command = list_rbind(list(!!!.x), names_to = ".id") |>
      zutils::separate_wider_dsv(
        .id,
        names(config_num_subjs),
        prefix = "patterns_group_stability"
      )
  ),
  tar_target(
    patterns_group_whole,
    arrow::read_parquet(file_cca_y) |>
      filter(time_id >= index_onset) |>
      calc_group_pattern()
  ),
  tar_target(
    data_igs_whole,
    calc_igs(patterns_indiv_whole, patterns_group_whole)
  ),
  tar_target(
    data_igs_partial_whole,
    calc_igs(
      patterns_indiv_whole |>
        mutate(pattern = map(pattern, get_resid, pattern_semantics)),
      patterns_group_whole |>
        mutate(pattern = map(pattern, get_resid, pattern_semantics))
    )
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
    # `na.rm` not supported in `open_dataset()`
    # https://github.com/apache/arrow/issues/44089
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
