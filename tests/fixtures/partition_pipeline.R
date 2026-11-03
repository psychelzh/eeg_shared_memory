library(targets)
tar_source()
tar_option_set(
  packages = "tidyverse",
  format = "qs",
  retrieval = "worker",
  resources = tar_resources(crew = tar_resources_crew(controller = "light")),
  controller = crew::crew_controller_group(
    crew::crew_controller_local(name = "light", workers = 2),
    crew::crew_controller_local(name = "heavy", workers = 1)
  )
)
heavy <- tar_resources(crew = tar_resources_crew(controller = "heavy"))
partition_sources <- c(
  patterns = "partitions_patterns",
  group = "partitions_group",
  isps_whole = "partitions_isps_whole",
  isps_dynamic = "partitions_isps_dynamic"
)
indices <- tibble::tibble(index_name = rlang::syms(c("dprime", "dprime_rem")))
shapes <- tibble::tibble(reference = rlang::syms(c("yd", "z")))
list(
  tar_target(patterns, {
    set.seed(47)
    crossing(cca_id = 1:3, subj_id = 1:8, time_id = seq(51L, 66L, 5L)) |>
      mutate(pattern = map(subj_id, ~ dist(matrix(rnorm(144), 12))))
  }),
  tar_target(group, {
    set.seed(48)
    mutate(patterns, pattern = map(pattern, ~ dist(matrix(rnorm(144), 12))))
  }),
  tar_target(yd, {
    set.seed(49)
    dist(matrix(rnorm(144), 12))
  }),
  tar_target(z, {
    set.seed(50)
    dist(matrix(rnorm(144), 12))
  }),
  tar_target(
    memory,
    tibble(subj_id = 1:8, dprime = seq_len(8), dprime_rem = -seq_len(8))
  ),
  tar_target(smc, {
    set.seed(51)
    dist(matrix(rnorm(64), 8))
  }),
  tar_target(z_smc, {
    set.seed(52)
    dist(matrix(rnorm(64), 8))
  }),
  tar_target(
    channel_regions,
    setNames(
      lapply(1:6, function(i) {
        as.integer(seq_len(12) %in% c(2 * i - 1, 2 * i))
      }),
      paste0("region", 1:6)
    )
  ),
  tar_target(
    raw_eeg,
    {
      set.seed(53)
      file <- file.path(targets::tar_path_store(), "test-eeg.qs2")
      qs2::qs_save(
        array(rnorm(12 * 61 * 8 * 6), c(12, 61, 8, 6)),
        file,
        nthreads = 1L
      )
      file
    },
    format = "file"
  ),
  tar_partition_registry(partition_sources, resources = heavy),
  tar_cluster_permutation(
    "iss_dynamic",
    data_expr = corr_patterns_1(patterns, yd, name = "iss"),
    data_perm_expr = corr_patterns_1(patterns, permute_dist(yd), name = "iss"),
    stats_expr = calc_stats_t(!!.x, iss),
    stats_perm_expr = calc_stats_t(!!.x, iss, alternative = "greater"),
    clusters_stats_expr = calc_clusters_stats(
      mutate(!!.x, p.value = convert_p2_p1(p.value, statistic)),
      !!.y
    ),
    reps = 2L,
    batches = 2L
  ),
  tar_cluster_permutation(
    "igs_dynamic",
    data_expr = corr_patterns(patterns, group, name = "igs"),
    data_perm_expr = corr_patterns(
      patterns,
      mutate(group, pattern = map(pattern, permute_dist)),
      name = "igs"
    ),
    stats_expr = calc_stats_t(!!.x, igs),
    stats_perm_expr = calc_stats_t(!!.x, igs, alternative = "greater"),
    reps = 2L,
    batches = 2L
  ),
  tar_cluster_permutation(
    "mantel_dynamic",
    data_expr = calc_mantel(
      filter(group, subj_id == 1L),
      yd,
      permutations = 19L
    ),
    data_perm_expr = calc_mantel(
      filter(group, subj_id == 1L),
      permute_dist(yd),
      permutations = 19L
    ),
    stats_expr = extract_stats_mantel(!!.x),
    stats_perm_expr = extract_stats_mantel(!!.x),
    reps = 2L,
    batches = 2L
  ),
  tarchetypes::tar_map(
    shapes,
    tar_cluster_permutation(
      "iws_dynamic",
      data_expr = corr_patterns_1(patterns, reference, name = "iws"),
      data_perm_expr = corr_patterns_1(
        patterns,
        permute_dist(reference),
        name = "iws"
      ),
      stats_expr = calc_stats_t(!!.x, iws),
      stats_perm_expr = calc_stats_t(!!.x, iws, alternative = "greater"),
      reps = 2L,
      batches = 2L
    )
  ),
  tar_partition_target(
    "isps_whole",
    filter(patterns, time_id == 51L) |> select(-time_id) |> calc_isps(),
    stream = FALSE
  ),
  tar_partition_target("isps_dynamic", calc_isps(patterns)),
  tar_partition_target(
    "isps_partial",
    regress_patterns_1(patterns, yd) |> calc_isps()
  ),
  tar_partition_replicates(
    "summary_isps_dynamic_permuted",
    patterns |>
      mutate(pattern = map(pattern, permute_dist)) |>
      calc_isps() |>
      summarise_isps(),
    reps = 4L,
    batch_size = 2L
  ),
  tar_mantel(
    "control",
    isps_whole,
    isps_dynamic,
    smc,
    z_smc,
    reps = 1L,
    batches = 2L,
    permutations = 19L
  ),
  tarchetypes::tar_map(
    indices,
    tar_cluster_permutation(
      "memory",
      corr_mem(data_iss_dynamic, memory, col_perf = index_name),
      corr_mem(
        data_iss_dynamic,
        mutate(memory, subj_id = sample(subj_id)),
        col_perf = index_name
      ),
      reps = 2L,
      batches = 2L
    )
  ),
  tar_target(
    partitions_eeg,
    write_eeg_region_partitions(
      raw_eeg,
      channel_regions,
      file.path(targets::tar_path_store(), "user", "eeg-regions")
    ),
    format = "file",
    resources = heavy
  ),
  tar_partition_target(
    "patterns_regions",
    calc_indiv_pattern_dynamic_region(.inputs$eeg),
    sources = c(eeg = "partitions_eeg")
  ),
  tar_partition_replicates(
    "probe",
    {
      start <- as.numeric(Sys.time())
      Sys.sleep(0.2)
      tibble(
        cca_id = unique(patterns$cca_id),
        time_id = unique(patterns$time_id),
        pid = Sys.getpid(),
        start = start,
        end = as.numeric(Sys.time())
      )
    },
    reps = 2L,
    batch_size = 1L
  )
)
