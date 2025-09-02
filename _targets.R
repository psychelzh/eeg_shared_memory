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
  garbage_collection = TRUE
)

future::plan(future.callr::callr)

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
    events_retrieval,
    "data/behav/retrieval.tsv",
    read = clean_events(!!.x, subjs)
  ),
  tar_target(mem_perf, calc_mem_perf(events_retrieval)),
  tar_target(mem_perf_precise, calc_mem_perf_precise(events_retrieval)),
  tarchetypes::tar_file_read(
    smc,
    "data/behav/simil.rds", # use pre-calculated
    read = readRDS(!!.x)$mat[[4]]
  ),
  tar_target(simil_mem, calc_simil_mem(mem_perf)),
  tar_target(memorability, calc_memorability(events_retrieval)),
  tar_target(
    memorability_content,
    calc_mem_content(events_retrieval, memorability)
  ),

  # representations (patterns) calculation ----
  ## configuration common to stimuli patterns ----
  tar_target(file_seq, "config/sem_sequence.mat", format = "file"),
  tar_target(
    mapping_word_trial,
    R.matlab::readMat(file_seq)$SM[, 1:2] |>
      as_tibble(.name_repair = ~ c("trial_id", "word_id"))
  ),
  tar_target(
    file_cca_y,
    "data/CorCAExtra/cca_y_subjs206.parquet",
    format = "file"
  ),
  # used for dynamic branches to save memory
  tar_target(subj_id_loop, seq_len(num_subj)),

  ## semantic patterns ----
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

  ## word shape patterns ----
  tar_target(
    file_charsim,
    "data/stimuli/word_shape_sims/char_similar.tsv",
    format = "file"
  ),
  tar_target(
    file_alexnet,
    "data/stimuli/word_shape_sims/alexnet.tsv",
    format = "file"
  ),
  tar_target(
    patterns_shapes,
    read_alexnet(file_alexnet, mapping_word_trial) |>
      add_row(
        model = "char_similar",
        layer = NA_character_,
        pattern = list(read_charsim(file_charsim, mapping_word_trial))
      )
  ),

  ## individualized patterns ----
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

  ## individualized patterns (orginal regions) ----
  tarchetypes::tar_file_read(
    channel_regions,
    "config/eeg_channel_labels_64.csv",
    read = read_csv(!!.x, show_col_types = FALSE)
  ),
  tar_target(
    file_eeg_data,
    "data-raw/grp_subjs206_nodemean_1000ms.qs2",
  ),
  tarchetypes::tar_map(
    list(region = 1:6),
    tar_target(
      patterns_indiv_dynamic_region,
      read_eeg_region(file_eeg_data, channel_regions, region) |>
        calc_indiv_pattern_dynamic_region()
    ),
    tar_target(
      patterns_indiv_whole_region,
      read_eeg_region(file_eeg_data, channel_regions, region) |>
        calc_indiv_pattern_region()
    ),
    tar_cluster_permutation(
      "iss_dynamic_region",
      data_expr = calc_iss(patterns_indiv_dynamic_region, pattern_semantics),
      data_perm_expr = calc_iss(
        patterns_indiv_dynamic_region,
        permute_dist(pattern_semantics)
      ),
      stats_expr = calc_stats_t(!!.x, iss, .by = time_id),
      stats_perm_expr = calc_stats_t(
        !!.x,
        iss,
        .by = time_id,
        alternative = "greater"
      ),
      clusters_stats_expr = calc_clusters_stats(
        mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
        !!.y,
        by = NULL
      )
    ),
    tar_target(
      data_iss_whole_region,
      calc_iss(patterns_indiv_whole_region, pattern_semantics)
    ),
    tar_target(
      stats_iss_whole_region,
      calc_stats_t(data_iss_whole_region, iss, .by = NULL)
    ),
    tar_cluster_permutation(
      "iws_dynamic_region",
      data_expr = calc_iws(patterns_indiv_dynamic_region, patterns_shapes),
      data_perm_expr = calc_iws(
        patterns_indiv_dynamic_region,
        patterns_shapes |>
          mutate(pattern = map(pattern, permute_dist))
      ),
      stats_expr = calc_stats_t(!!.x, iws, .by = c(time_id, model, layer)),
      stats_perm_expr = calc_stats_t(
        !!.x,
        iws,
        .by = c(time_id, model, layer),
        alternative = "greater"
      ),
      clusters_stats_expr = calc_clusters_stats(
        mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
        !!.y,
        by = c("model", "layer")
      )
    ),
    tar_target(
      data_iws_whole_region,
      calc_iws(patterns_indiv_whole_region, patterns_shapes)
    ),
    tar_target(
      stats_iws_whole_region,
      calc_stats_t(data_iws_whole_region, iws, .by = c(model, layer))
    )
  ),

  ## group averaged patterns ----
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

  # IGS: individual to group averaged pattern similarity ----
  ## IGS calculation ----
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
  tar_target(stats_igs_whole, calc_stats_t(data_igs_whole, igs, .by = cca_id)),
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
  tar_cluster_permutation(
    "igs_dynamic",
    data_expr = calc_igs(patterns_indiv_dynamic, patterns_group_dynamic_loo),
    data_perm_expr = calc_igs(
      patterns_indiv_dynamic,
      patterns_group_dynamic_loo |>
        mutate(pattern = map(pattern, permute_dist))
    ),
    stats_expr = calc_stats_t(!!.x, igs),
    stats_perm_expr = calc_stats_t(!!.x, igs, alternative = "greater"),
    clusters_stats_expr = calc_clusters_stats(
      mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
      !!.y
    )
  ),
  tar_target(igs_comparison, compare_inter_patterns(data_igs_whole)),

  ## IGS predicts memory ----
  tar_target(stats_igs_mem_whole, correlate_mem_perf(data_igs_whole, mem_perf)),
  tar_cluster_permutation(
    "igs_mem_dynamic",
    correlate_mem_perf(data_igs_dynamic, mem_perf),
    correlate_mem_perf(
      data_igs_dynamic,
      mutate(mem_perf, subj_id = sample(subj_id)),
      alternative = "greater"
    ),
    clusters_stats_expr = calc_clusters_stats(
      mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
      !!.y
    )
  ),
  # more precise memory performance
  tarchetypes::tar_map(
    config_mem_precise,
    tar_target(
      stats_igs_mem_whole,
      correlate_mem_perf(
        data_igs_whole,
        select(mem_perf_precise, subj_id, dprime = index_name)
      )
    ),
    tar_cluster_permutation(
      "igs_mem_dynamic",
      correlate_mem_perf(
        data_igs_dynamic,
        select(mem_perf_precise, subj_id, dprime = index_name)
      ),
      correlate_mem_perf(
        data_igs_dynamic,
        mem_perf_precise |>
          mutate(subj_id = sample(subj_id)) |>
          select(subj_id, dprime = index_name),
        alternative = "greater"
      ),
      clusters_stats_expr = calc_clusters_stats(
        mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
        !!.y
      )
    )
  ),

  # GSS: group averaged patterns and semantic pattern ----
  tar_target(
    data_gss_whole,
    calc_mantel(patterns_group_whole, pattern_semantics)
  ),
  tar_target(stats_gss_whole, extract_stats_mantel(data_gss_whole)),
  tar_cluster_permutation(
    "gss_dynamic",
    data_expr = calc_mantel(patterns_group_dynamic, pattern_semantics),
    data_perm_expr = calc_mantel(
      patterns_group_dynamic,
      permute_dist(pattern_semantics)
    ),
    stats_expr = extract_stats_mantel(!!.x),
    stats_perm_expr = extract_stats_mantel(!!.x)
  ),

  # GWS: group averaged and word shape ----
  tar_target(
    data_gws_whole,
    calc_mantel2(patterns_group_whole, patterns_shapes)
  ),
  tar_target(stats_gws_whole, extract_stats_mantel(data_gws_whole)),
  tar_cluster_permutation(
    "gws_dynamic",
    data_expr = calc_mantel2(patterns_group_dynamic, patterns_shapes),
    data_perm_expr = calc_mantel2(
      patterns_group_dynamic,
      patterns_shapes |>
        mutate(pattern = map(pattern, permute_dist))
    ),
    stats_expr = extract_stats_mantel(!!.x),
    stats_perm_expr = extract_stats_mantel(!!.x),
    clusters_stats_expr = calc_clusters_stats(
      !!.x,
      !!.y,
      by = c("cca_id", "model", "layer")
    ),
    reps = 10, # to save memory
    batches = 100
  ),

  # ISS: individual to semantic patterns ----
  ## ISS calculation ----
  tar_cluster_permutation(
    "iss_dynamic",
    data_expr = calc_iss(patterns_indiv_dynamic, pattern_semantics),
    data_perm_expr = calc_iss(
      patterns_indiv_dynamic,
      permute_dist(pattern_semantics)
    ),
    stats_expr = calc_stats_t(!!.x, iss),
    stats_perm_expr = calc_stats_t(!!.x, iss, alternative = "greater"),
    clusters_stats_expr = calc_clusters_stats(
      mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
      !!.y
    )
  ),
  tar_target(data_iss_whole, calc_iss(patterns_indiv_whole, pattern_semantics)),
  tar_target(stats_iss_whole, calc_stats_t(data_iss_whole, iss, .by = cca_id)),
  tar_target(iss_comparison, compare_inter_patterns(data_iss_whole)),

  ## ISS predicts memory ----
  tar_target(stats_iss_mem_whole, correlate_mem_perf(data_iss_whole, mem_perf)),
  tar_target(comparison_iss_mem, compare_iss_mem(stats_iss_mem_whole)),
  tar_cluster_permutation(
    "iss_mem_dynamic",
    correlate_mem_perf(data_iss_dynamic, mem_perf),
    correlate_mem_perf(
      data_iss_dynamic,
      mutate(mem_perf, subj_id = sample(subj_id)),
      alternative = "greater"
    ),
    clusters_stats_expr = calc_clusters_stats(
      mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
      !!.y
    )
  ),
  tar_target(iss_mem_combine, fit_mem_pred(mem_perf, data_iss_whole)),
  tar_target(
    iss_r2_combine,
    calc_iss_ind_r2(patterns_indiv_whole, pattern_semantics)
  ),
  tar_target(
    iss_r2_mem_combine,
    caret::train(
      dprime ~ r2,
      data = inner_join(iss_r2_combine, mem_perf, by = "subj_id"),
      method = "lm",
      trControl = caret::trainControl(method = "LOOCV")
    )
  ),
  # more precise memory performance
  tarchetypes::tar_map(
    config_mem_precise,
    tar_target(
      stats_iss_mem_whole,
      correlate_mem_perf(
        data_iss_whole,
        select(mem_perf_precise, subj_id, dprime = index_name)
      )
    ),
    tar_cluster_permutation(
      "iss_mem_dynamic",
      correlate_mem_perf(
        data_iss_dynamic,
        select(mem_perf_precise, subj_id, dprime = index_name)
      ),
      correlate_mem_perf(
        data_iss_dynamic,
        mem_perf_precise |>
          mutate(subj_id = sample(subj_id)) |>
          select(subj_id, dprime = index_name),
        alternative = "greater"
      ),
      clusters_stats_expr = calc_clusters_stats(
        mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
        !!.y
      )
    )
  ),

  # mediation analysis ----
  tar_target(
    data_med,
    data_igs_whole |>
      inner_join(data_iss_whole, by = c("subj_id", "cca_id")) |>
      inner_join(mem_perf, by = "subj_id")
  ),
  tar_target(
    cor_igs_iss_whole,
    data_med |>
      summarise(broom::tidy(cor.test(iss, igs)), .by = cca_id)
  ),
  tarchetypes::tar_file_read(
    model_med,
    "config/mediation.lav",
    read = readr::read_lines(!!.x)
  ),
  tar_target(data_combined, combine_data_ccas(data_med)),
  tar_target(
    fit_med_iss_igs_dprime,
    fit_med(model_med, data_combined, X = "iss", Y = "dprime", M = "igs")
  ),
  tar_target(
    fit_med_igs_iss_dprime,
    fit_med(model_med, data_combined, X = "igs", Y = "dprime", M = "iss")
  ),

  # control analyses ----
  ## IWS: individual patterns and word shape pattern ----
  ### IWS calculation ----
  tar_cluster_permutation(
    "iws_dynamic",
    data_expr = calc_iws(patterns_indiv_dynamic, patterns_shapes),
    data_perm_expr = calc_iws(
      patterns_indiv_dynamic,
      patterns_shapes |>
        mutate(pattern = map(pattern, permute_dist))
    ),
    stats_expr = calc_stats_t(
      !!.x,
      iws,
      .by = c(cca_id, time_id, model, layer)
    ),
    stats_perm_expr = calc_stats_t(
      !!.x,
      iws,
      .by = c(cca_id, time_id, model, layer),
      alternative = "greater"
    ),
    clusters_stats_expr = calc_clusters_stats(
      mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
      !!.y,
      by = c("cca_id", "model", "layer")
    ),
    reps = 100, # to save memory
    batches = 10
  ),
  tar_target(data_iws_whole, calc_iws(patterns_indiv_whole, patterns_shapes)),
  tar_target(
    stats_iws_whole,
    calc_stats_t(data_iws_whole, iws, .by = c(cca_id, model, layer))
  ),
  tar_target(
    iws_comparison,
    data_iws_whole |>
      nest(.by = c(model, layer)) |>
      mutate(
        comparison = map(data, compare_inter_patterns),
        .keep = "unused"
      ) |>
      unnest(comparison)
  ),

  ### IWS predicts memory ----
  tar_target(stats_iws_mem_whole, correlate_mem_perf(data_iws_whole, mem_perf)),
  tar_target(comparison_iws_mem, compare_iss_mem(stats_iws_mem_whole)),
  tar_cluster_permutation(
    "iws_mem_dynamic",
    correlate_mem_perf(data_iws_dynamic, mem_perf),
    correlate_mem_perf(
      data_iws_dynamic,
      mutate(mem_perf, subj_id = sample(subj_id))
    ),
    clusters_stats_expr = calc_clusters_stats(
      mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
      mutate(!!.y, p.value = convert_p2_p1(p.value, statistic)),
      by = c("cca_id", "model", "layer")
    )
  ),
  tar_target(
    clusters_stats_less_iws_mem_dynamic,
    calc_clusters_stats(
      stats_iws_mem_dynamic |>
        mutate(
          p.value = convert_p2_p1(
            p.value,
            statistic,
            alternative = "less"
          )
        ),
      stats_iws_mem_dynamic_permuted |>
        mutate(
          p.value = convert_p2_p1(
            p.value,
            statistic,
            alternative = "less"
          )
        ),
      alternative = "less",
      by = c("cca_id", "model", "layer")
    )
  ),

  ## regress semantic from group averaged ----
  tar_target(
    data_igs_partial_whole,
    calc_igs(
      patterns_indiv_whole |>
        mutate(pattern = map(pattern, get_resid, pattern_semantics)),
      patterns_group_whole_loo |>
        mutate(pattern = map(pattern, get_resid, pattern_semantics))
    )
  ),
  tar_target(
    igs_comp_partial,
    compare_partial(data_igs_whole, data_igs_partial_whole)
  ),
  tar_target(
    lm_mem_igs_partial,
    fit_mem_pred(mem_perf, data_igs_partial_whole)
  ),
  tar_target(
    lm_mem_iss_igs_partial,
    fit_mem_pred(mem_perf, data_igs_partial_whole, data_iss_whole)
  ),
  tar_target(
    data_igs_partial_dynamic,
    calc_igs(
      patterns_indiv_dynamic |>
        mutate(pattern = map(pattern, get_resid, pattern_semantics)),
      patterns_group_dynamic_loo |>
        mutate(pattern = map(pattern, get_resid, pattern_semantics))
    )
  ),
  tar_cluster_permutation(
    "igs_partial_mem_dynamic",
    correlate_mem_perf(data_igs_partial_dynamic, mem_perf),
    correlate_mem_perf(
      data_igs_partial_dynamic,
      mutate(mem_perf, subj_id = sample(subj_id)),
      alternative = "greater"
    ),
    clusters_stats_expr = calc_clusters_stats(
      mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
      !!.y
    )
  ),

  ## regress group averaged from semantic patterns ----
  tar_target(
    data_iss_partial_whole,
    calc_iss2(
      patterns_indiv_whole |>
        inner_join(
          patterns_group_whole_loo,
          by = c("subj_id", "cca_id"),
          suffix = c("_indiv", "_group")
        ) |>
        mutate(
          pattern = map2(pattern_indiv, pattern_group, get_resid),
          .keep = "unused"
        ),
      patterns_group_whole_loo |>
        mutate(pattern = map(pattern, \(x) get_resid(pattern_semantics, x)))
    )
  ),
  tar_target(
    iss_comp_partial,
    compare_partial(data_iss_whole, data_iss_partial_whole)
  ),
  tar_target(
    lm_mem_iss_partial,
    fit_mem_pred(mem_perf, data_iss_partial_whole)
  ),
  tar_target(
    lm_mem_igs_iss_partial,
    fit_mem_pred(mem_perf, data_iss_partial_whole, data_igs_whole)
  ),
  tar_target(
    data_iss_partial_dynamic,
    calc_iss2(
      patterns_indiv_dynamic |>
        inner_join(
          patterns_group_dynamic_loo,
          by = c("subj_id", "cca_id", "time_id"),
          suffix = c("_indiv", "_group")
        ) |>
        mutate(
          pattern = map2(pattern_indiv, pattern_group, get_resid),
          .keep = "unused"
        ),
      patterns_group_dynamic_loo |>
        mutate(pattern = map(pattern, \(x) get_resid(pattern_semantics, x)))
    )
  ),
  tar_cluster_permutation(
    "iss_partial_mem_dynamic",
    correlate_mem_perf(data_iss_partial_dynamic, mem_perf),
    correlate_mem_perf(
      data_iss_partial_dynamic,
      mutate(mem_perf, subj_id = sample(subj_id)),
      alternative = "greater"
    ),
    clusters_stats_expr = calc_clusters_stats(
      mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
      !!.y
    )
  ),

  # intersubject pattern similarity (ISPS) ----
  tar_target(data_isps_whole, calc_isps(patterns_indiv_whole)),
  tarchetypes::tar_rep(
    data_isps_whole_permuted,
    patterns_indiv_whole |>
      mutate(pattern = map(pattern, permute_dist)) |>
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
  tar_target(
    summary_isps_dynamic,
    summarise_isps(data_isps_dynamic, se = TRUE)
  ),
  tarchetypes::tar_rep(
    data_isps_dynamic_permuted,
    patterns_indiv_dynamic |>
      mutate(pattern = map(pattern, permute_dist)) |>
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
    data_perm_expr = calc_mantel(data_isps_dynamic, permute_dist(smc)),
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
    "sync_smc_dynamic",
    data_expr = calc_mantel(sync_dynamic, smc),
    data_perm_expr = calc_mantel(sync_dynamic, permute_dist(smc)),
    stats_expr = extract_stats_mantel(!!.x),
    stats_perm_expr = extract_stats_mantel(!!.x)
  ),

  # WIP: compare ISPS and synchronization in predicting SMC ----
  tar_target(
    mantel_isps_sync,
    data_isps_whole |>
      inner_join(sync_whole, by = "cca_id") |>
      mutate(
        mantel = map2(
          isps,
          pattern,
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
          isps,
          pattern,
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
