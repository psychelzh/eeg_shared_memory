library(targets)

tar_option_set(
  packages = c("tidyverse"),
  controller = if (Sys.info()["nodename"] %in% c("shadow", "hippocampus")) {
    crew.cluster::crew_controller_sge(
      name = "sge",
      workers = 25,
      seconds_idle = 30
    )
  } else {
    crew::crew_controller_local(
      name = "local",
      workers = 12
    )
  },
  garbage_collection = TRUE,
  memory = "transient",
  retrieval = "worker",
  storage = "worker"
)

tar_source()

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
        arrow::open_dataset(file_cca_y) |>
          filter(time_id >= index_onset, subj_id %in% subjs) |>
          collect() |>
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
  ),
  tarchetypes::tar_rep2(
    data_gss_whole_resampled,
    patterns_group_whole_resampled |>
      mutate(
        gss = map_dbl(pattern, cor, pattern_semantics),
        .keep = "unused"
      ),
    patterns_group_whole_resampled
  )
)

list(
  # behavioral data ----
  tarchetypes::tar_file_read(
    subjs,
    "data/subj_206.txt",
    read = scan(!!.x)
  ),
  tarchetypes::tar_file_read(
    mem_perf,
    "data/behav/retrieval.tsv",
    read = calc_mem_perf(read_tsv(!!.x, show_col_types = FALSE), subjs)
  ),
  tarchetypes::tar_file_read(
    smc,
    "data/behav/simil.rds", # use pre-calculated
    read = readRDS(!!.x)$mat[[4]]
  ),
  tar_target(simil_mem, calc_simil_mem(mem_perf)),

  # stimuli patterns ----
  tar_target(file_seq, "config/sem_sequence.mat", format = "file"),
  tar_target(
    mapping_word_trial,
    raveio::read_mat(file_seq)$SM[, 1:2] |>
      as_tibble(.name_repair = ~ c("trial_id", "word_id"))
  ),
  tar_target(file_w2v, "data/stimuli/words_w2v.txt", format = "file"),
  tar_target(
    pattern_semantics,
    mapping_word_trial |>
      inner_join(
        read_table(file_w2v, show_col_types = FALSE, col_names = FALSE),
        by = c("word_id" = "X1")
      ) |>
      filter(trial_id > 0) |>
      select(-word_id, -X2) |>
      column_to_rownames("trial_id") |>
      proxy::simil(method = "cosine")
  ),
  tar_target(
    file_word_shape,
    "data/stimuli/words_shape_similarity.tsv",
    format = "file"
  ),
  tar_target(
    pattern_shapes,
    {
      order <- with(mapping_word_trial, word_id[trial_id > 0])
      x <- read_tsv(file_word_shape, show_col_types = FALSE)$similarity |>
        pracma::squareform()
      as.dist(x[order, order])
    }
  ),
  tar_target(
    file_cca_y,
    "data/CorCAExtra/cca_y_subjs206.parquet",
    format = "file"
  ),
  tar_target(subj_id_loop, seq_len(num_subj)),

  # individualized patterns ----
  tar_target(
    patterns_indiv_dynamic,
    arrow::open_dataset(file_cca_y) |>
      filter(subj_id == subj_id_loop) |>
      collect() |>
      calc_indiv_pattern_dynamic(),
    pattern = map(subj_id_loop)
  ),
  tar_target(
    patterns_indiv_whole,
    arrow::open_dataset(file_cca_y) |>
      filter(time_id >= index_onset) |>
      collect() |>
      calc_indiv_pattern()
  ),

  # group averaged patterns ----
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
  tarchetypes::tar_combine(
    data_gss_whole_resampled,
    targets_patterns_group_whole_resampled$data_gss_whole_resampled,
    command = list_rbind(list(!!!.x), names_to = ".id") |>
      zutils::separate_wider_dsv(
        .id,
        names(config_num_subjs),
        prefix = "data_gss_whole_resampled"
      )
  ),
  tar_target(
    patterns_group_dynamic,
    arrow::read_parquet(file_cca_y) |>
      calc_group_pattern_dynamic()
  ),
  tar_target(
    patterns_group_whole,
    arrow::open_dataset(file_cca_y) |>
      filter(time_id >= index_onset) |>
      collect() |>
      calc_group_pattern()
  ),

  # individual to group averaged pattern similarity ----
  tar_target(
    # leave one out
    patterns_group_whole_loo,
    arrow::open_dataset(file_cca_y) |>
      filter(time_id >= index_onset, subj_id != subj_id_loop) |>
      collect() |>
      calc_group_pattern() |>
      add_column(subj_id = subj_id_loop, .before = 1L),
    pattern = map(subj_id_loop)
  ),
  tar_target(
    data_igs_whole,
    calc_igs(patterns_indiv_whole, patterns_group_whole_loo)
  ),
  tar_target(
    # leave one out
    patterns_group_dynamic_loo,
    arrow::open_dataset(file_cca_y) |>
      filter(subj_id != subj_id_loop) |>
      collect() |>
      calc_group_pattern_dynamic() |>
      add_column(subj_id = subj_id_loop, .before = 1L),
    pattern = map(subj_id_loop)
  ),
  tar_target(
    data_igs_dynamic,
    calc_igs(patterns_indiv_dynamic, patterns_group_dynamic_loo)
  ),

  # IGS predicts memory ----
  tar_target(stats_igs_mem_whole, calc_igs_mem(data_igs_whole, mem_perf)),
  tar_cluster_permutation(
    "igs_mem_dynamic",
    calc_igs_mem(data_igs_dynamic, mem_perf),
    calc_igs_mem(
      data_igs_dynamic,
      mutate(mem_perf, subj_id = sample(subj_id)),
      alternative = "greater"
    ),
    clusters_stats_expr = calc_clusters_stats(
      mutate(!!.x, p.value = convert_p2_p1(statistic, p.value)),
      !!.y
    )
  ),

  # group averaged patterns and semantic pattern ----
  tar_target(data_gss_whole, calc_mantel(patterns_group_whole, pattern_semantics)),
  tar_target(stats_gss_whole, extract_stats_mantel(data_gss_whole)),
  tar_cluster_permutation(
    "gss_dynamic",
    data_expr = calc_mantel(patterns_group_dynamic, pattern_semantics),
    data_perm_expr = calc_mantel(
      patterns_group_dynamic,
      seriation::permute(pattern_semantics, sample.int(150L))
    ),
    stats_expr = extract_stats_mantel(!!.x),
    stats_perm_expr = extract_stats_mantel(!!.x)
  ),
  # gcs: group averaged and character similarity
  tar_target(data_gcs_whole, calc_mantel(patterns_group_whole, pattern_shapes)),
  tar_target(stats_gcs_whole, extract_stats_mantel(data_gcs_whole)),

  # individual patterns and semantic pattern similarity (ISS) ----
  tar_cluster_permutation(
    "iss_dynamic",
    data_expr = calc_iss(patterns_indiv_dynamic, pattern_semantics),
    data_perm_expr = calc_iss(
      patterns_indiv_dynamic,
      seriation::permute(pattern_semantics, sample.int(150L))
    ),
    stats_expr = calc_iss_stats(!!.x),
    stats_perm_expr = calc_iss_stats(!!.x, alternative = "greater"),
    clusters_stats_expr = calc_clusters_stats(
      mutate(!!.x, p.value = convert_p2_p1(statistic, p.value)),
      !!.y
    )
  ),
  tar_target(data_iss_whole, calc_iss(patterns_indiv_whole, pattern_semantics)),
  tar_target(stats_iss_whole, calc_iss_stats(data_iss_whole, .by = cca_id)),
  tar_target(iss_comparison, compare_iss(data_iss_whole)),

  # ISS predicts memory ----
  tar_target(stats_iss_mem_whole, calc_iss_mem(data_iss_whole, mem_perf)),
  tar_target(comparison_iss_mem, compare_iss_mem(stats_iss_mem_whole)),
  tar_cluster_permutation(
    "iss_mem_dynamic",
    calc_iss_mem(data_iss_dynamic, mem_perf),
    calc_iss_mem(
      data_iss_dynamic,
      mutate(mem_perf, subj_id = sample(subj_id)),
      alternative = "greater"
    ),
    clusters_stats_expr = calc_clusters_stats(
      mutate(!!.x, p.value = convert_p2_p1(statistic, p.value)),
      !!.y
    )
  ),

  # regress semantic from group averaged ----
  tar_target(
    data_igs_partial_whole,
    calc_igs(
      patterns_indiv_whole |>
        mutate(pattern = map(pattern, get_resid, pattern_semantics)),
      patterns_group_whole_loo |>
        mutate(pattern = map(pattern, get_resid, pattern_semantics))
    )
  ),
  tar_target(igs_comp_partial, compare_igs(data_igs_whole, data_igs_partial_whole)),
  tar_target(lm_mem_igs_partial, fit_mem_pred(mem_perf, data_igs_partial_whole)),
  tar_target(
    lm_mem_iss_igs_partial,
    fit_mem_pred(mem_perf, data_igs_partial_whole, data_iss_whole)
  ),

  # intersubject pattern similarity (ISPS) ----
  tar_target(data_isps_whole, calc_isps(patterns_indiv_whole)),
  tarchetypes::tar_rep(
    data_isps_whole_permuted,
    patterns_indiv_whole |>
      mutate(
        pattern = map(
          pattern,
          \(x) seriation::permute(x, sample.int(150L))
        )
      ) |>
      calc_isps(),
    reps = 10,
    batches = 100
  ),
  tarchetypes::tar_rep2(
    summary_isps_whole_permuted,
    summarise_isps(data_isps_whole_permuted),
    data_isps_whole_permuted
  ),
  tar_target(
    stats_isps_whole,
    summarise_isps(data_isps_whole, se = TRUE) |>
      calc_stats_isps(summary_isps_whole_permuted)
  ),
  tar_target(data_isps_dynamic, calc_isps(patterns_indiv_dynamic)),
  tar_target(summary_isps_dynamic, summarise_isps(data_isps_dynamic, se = TRUE)),
  tarchetypes::tar_rep(
    data_isps_dynamic_permuted,
    patterns_indiv_dynamic |>
      mutate(
        pattern = map(
          pattern,
          \(x) seriation::permute(x, sample.int(150L))
        )
      ) |>
      calc_isps(),
    reps = 10,
    batches = 100
  ),
  tarchetypes::tar_rep2(
    summary_isps_dynamic_permuted,
    summarise_isps(data_isps_dynamic_permuted),
    data_isps_dynamic_permuted
  ),
  tar_target(
    stats_isps_dynamic,
    calc_stats_isps(summary_isps_dynamic, summary_isps_dynamic_permuted)
  ),
  tar_target(
    stats_isps_dynamic_permuted,
    summary_isps_dynamic_permuted |>
      mutate(
        # TODO: find better method for p value calculations in future
        p_perm = percent_rank(desc(isps_mean)),
        .by = c(cca_id, time_id)
      )
  ),
  tar_target(
    clusters_stats_isps_dynamic,
    calc_clusters_stats(
      stats_isps_dynamic,
      stats_isps_dynamic_permuted,
      col_statistic = isps_mean,
      col_p_value = p_perm
    )
  ),

  # ISPS and shared memory content (SMC) ----
  tar_target(data_isps_smc_whole, calc_mantel(data_isps_whole, smc)),
  tar_target(stats_isps_smc_whole, extract_stats_mantel(data_isps_smc_whole)),
  tar_cluster_permutation(
    "isps_smc_dynamic",
    data_expr = calc_mantel(data_isps_dynamic, smc),
    data_perm_expr = calc_mantel(
      data_isps_dynamic,
      seriation::permute(smc, sample.int(206L))
    ),
    stats_expr = extract_stats_mantel(!!.x),
    stats_perm_expr = extract_stats_mantel(!!.x)
  ),
  # control for memory ability
  tar_mantel(
    "isps_smc_partial_ability",
    data_isps_whole,
    data_isps_dynamic,
    smc,
    simil_mem
  ),
  # control for group-averaged representation
  tar_target(
    data_isps_partial_group_whole,
    regress_patterns(patterns_indiv_whole, patterns_group_whole) |>
      calc_isps()
  ),
  tar_target(
    data_isps_partial_group_dynamic,
    regress_patterns(patterns_indiv_dynamic, patterns_group_dynamic) |>
      calc_isps()
  ),
  tar_mantel(
    "isps_smc_partial_group",
    data_isps_partial_group_whole,
    data_isps_partial_group_dynamic,
    smc
  ),
  tar_mantel(
    "isps_smc_partial_group_ability",
    data_isps_partial_group_whole,
    data_isps_partial_group_dynamic,
    smc,
    simil_mem
  ),
  # control for semantic representation
  tar_target(
    data_isps_partial_semantic_whole,
    patterns_indiv_whole |>
      mutate(pattern = map(pattern, regress_pattern, pattern_semantics)) |>
      calc_isps()
  ),
  tar_target(
    data_isps_partial_semantic_dynamic,
    patterns_indiv_dynamic |>
      mutate(pattern = map(pattern, regress_pattern, pattern_semantics)) |>
      calc_isps()
  ),
  tar_mantel(
    "isps_smc_partial_semantic",
    data_isps_partial_semantic_whole,
    data_isps_partial_semantic_dynamic,
    smc
  ),
  tar_mantel(
    "isps_smc_partial_semantic_ability",
    data_isps_partial_semantic_whole,
    data_isps_partial_semantic_dynamic,
    smc,
    simil_mem
  ),

  # shared and individualized patterns ----
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
    arrow::read_parquet(file_cca_y) |>
      mutate(half = if_else(trial_id <= 75, "first", "second")) |>
      summarise(
        y_avg = mean(y, na.rm = TRUE),
        .by = c(subj_id, cca_id, time_id, half)
      )
  ),
  tar_target(
    sync_within_halves,
    cca_y_halves |>
      filter(time_id >= index_onset) |>
      calc_sync_within_halves()
  ),
  tar_target(
    sync_between_halves,
    cca_y_halves |>
      filter(time_id >= index_onset) |>
      calc_sync_between_halves()
  ),
  tar_target(sync_inter_intra, prepare_sync_inter_intra(sync_between_halves)),

  # inter-subject synchronization predicts memory ----
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
  tar_target(
    sync_whole,
    whole_erps |>
      filter(time_id >= index_onset) |>
      calc_sync_whole()
  ),
  tar_target(sync_smc_whole, calc_mantel(sync_whole, smc)),
  tar_target(
    stats_sync_smc_whole,
    extract_stats_mantel(sync_smc_whole)
  ),
  tar_target(sync_dynamic, calc_sync_dynamic(whole_erps)),
  # the verbose names are for the compatibility with history relics
  tar_cluster_permutation(
    data_expr = calc_mantel(sync_dynamic, smc),
    data_perm_expr = calc_mantel(
      sync_dynamic,
      seriation::permute(smc, sample.int(206L))
    ),
    data_name = "sync_smc_dynamic",
    stats_expr = extract_stats_mantel(!!.x),
    stats_perm_expr = extract_stats_mantel(!!.x),
    stats_name = "stats_sync_smc_dynamic",
    clusters_stats_name = "clusters_stats_sync_smc_dynamic"
  ),

  # WIP: compare ISPS and synchronization in predicting SMC ----
  tar_target(
    mantel_isps_sync,
    data_isps_whole |>
      inner_join(sync_whole, by = "cca_id") |>
      mutate(
        mantel = map2(
          isps, pattern,
          vegan::mantel,
          permutations = 9999
        ),
        .keep = "unused"
      )
  ),
  tar_target(stats_mantel_isps_sync, extract_stats_mantel(mantel_isps_sync)),
  tar_target(
    fit_smc_isps_sync,
    data_isps_whole |>
      inner_join(sync_whole, by = "cca_id") |>
      mutate(
        fit = map2(
          isps, pattern,
          \(isps, sync) lm(smc ~ isps + sync)
        ),
        .keep = "unused"
      )
  ),
  tar_target(
    stats_smc_isps_sync,
    fit_smc_isps_sync |>
      reframe(
        map(fit, broom::tidy) |>
          list_rbind(),
        .by = cca_id
      )
  )
)
