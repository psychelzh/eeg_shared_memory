# Opt-in actual-worker integration, all generated inputs and a temporary store.
library(targets)
library(tidyverse)
tar_source()
script <- "tests/fixtures/partition_pipeline.R"
store <- tempfile("partition-pipeline-", tmpdir = Sys.getenv("TEMP"))
legacy <- tempfile(fileext = ".R")
lines <- readLines(script)
boundary <- grep("^  tar_partition_registry", lines)[1L]
writeLines(
  c(
    lines[seq_len(boundary - 1L)],
    'tar_target(data_iss_dynamic,corr_patterns_1(patterns,yd,name="iss")),',
    'tar_target(partitions_patterns,write_table_partitions(patterns,"cca_id",',
    ' file.path(targets::tar_path_store(),"user","legacy-patterns")),format="file"),',
    'tar_target(jobs_stats_iss_dynamic_permuted,',
    ' make_partition_jobs(list(patterns=partitions_patterns))),',
    'tar_target(plan_stats_iss_dynamic_permuted,',
    ' make_permutation_plan(4L,2L,targets::tar_seed_get())),',
    'tar_target(batches_stats_iss_dynamic_permuted,1:2),',
    'tar_target(stats_iss_dynamic_permuted,run_partition_batch(',
    ' jobs_stats_iss_dynamic_permuted,plan_stats_iss_dynamic_permuted,',
    ' batches_stats_iss_dynamic_permuted,',
    ' function(.inputs,.replicate_seed) calc_stats_t(',
    ' corr_patterns_1(.inputs$patterns,permute_dist(yd),name="iss"),',
    ' iss,alternative="greater")),',
    ' pattern=cross(jobs_stats_iss_dynamic_permuted,',
    ' batches_stats_iss_dynamic_permuted)))'
  ),
  legacy
)
tar_make(script = legacy, store = store, reporter = "silent")
legacy_files <- list.files(file.path(store, "user"), recursive = TRUE)
legacy_null <- tar_read_raw("stats_iss_dynamic_permuted", store = store)
tar_make(script = script, store = store, reporter = "balanced")
patterns <- tar_read_raw("patterns", store = store)
group <- tar_read_raw("group", store = store)
yd <- tar_read_raw("yd", store = store)
source_map <- c(patterns = "unused", group = "unused")
inputs <- list(patterns = patterns, group = group)
sort_result <- function(x) {
  x |> arrange(across(any_of(c("cca_id", "tar_rep", "time_id", "subj_id"))))
}
check_null <- function(name, expr) {
  plan <- tar_read_raw(paste0("plan_", name), store = store)
  fun <- eval(partition_callback(expr, source_map))
  expected <- map(seq_len(nrow(plan) - 1L) + 1L, function(i) {
    with_permutation_seed(plan$seed[i], fun(inputs, plan$seed[i])) |>
      mutate(tar_batch = plan$tar_batch[i], tar_rep = plan$tar_rep[i])
  }) |>
    list_rbind() |>
    sort_result()
  actual <- tar_read_raw(name, store = store) |> sort_result()
  stopifnot(isTRUE(all.equal(expected, actual, tolerance = 1e-12)))
}
check_null(
  "stats_iss_dynamic_permuted",
  quote(calc_stats_t(
    corr_patterns_1(patterns, permute_dist(yd), name = "iss"),
    iss,
    alternative = "greater"
  ))
)
check_null(
  "stats_igs_dynamic_permuted",
  quote(calc_stats_t(
    corr_patterns(
      patterns,
      mutate(group, pattern = map(pattern, permute_dist)),
      name = "igs"
    ),
    igs,
    alternative = "greater"
  ))
)
check_null(
  "stats_mantel_dynamic_permuted",
  quote(extract_stats_mantel(calc_mantel(
    filter(group, subj_id == 1L),
    permute_dist(yd),
    permutations = 19L
  )))
)
check_null(
  "summary_isps_dynamic_permuted",
  quote(
    patterns |>
      mutate(pattern = map(pattern, permute_dist)) |>
      calc_isps() |>
      summarise_isps()
  )
)
actual_isps <- tar_read_raw("isps_dynamic", store = store) |> sort_result()
stopifnot(isTRUE(all.equal(actual_isps, sort_result(calc_isps(patterns)))))
stopifnot(isTRUE(all.equal(
  sort_result(legacy_null),
  sort_result(tar_read_raw("stats_iss_dynamic_permuted", store = store)),
  tolerance = 1e-12
)))
memory_a <- tar_read_raw("stats_memory_dprime", store = store)
memory_b <- tar_read_raw("stats_memory_dprime_rem", store = store)
stopifnot(isTRUE(all.equal(memory_a$statistic, -memory_b$statistic)))
shape_a <- tar_read_raw("stats_iws_dynamic_yd", store = store)
shape_b <- tar_read_raw("stats_iws_dynamic_z", store = store)
stopifnot(
  !isTRUE(all.equal(shape_a$statistic, shape_b$statistic)),
  nrow(tar_read_raw("stats_iws_dynamic_permuted_yd", store = store)) ==
    4L * 3L * 4L
)
control <- tar_read_raw("stats_control_dynamic_permuted", store = store)
stopifnot(nrow(control) == 2L * 3L * 4L)
regions <- tar_read_raw("patterns_regions", store = store)
reference <- read_eeg_regions(
  tar_read_raw("raw_eeg", store = store),
  tar_read_raw("channel_regions", store = store)
) |>
  mutate(data = map(data, calc_indiv_pattern_dynamic_region)) |>
  unnest(data)
stopifnot(isTRUE(all.equal(regions, reference)))
probe <- tar_read_raw("probe", store = store) |> distinct(pid, start, end)
overlap <- outer(probe$start, probe$end, `<`) &
  t(outer(probe$start, probe$end, `<`)) &
  outer(probe$pid, probe$pid, `!=`)
stopifnot(n_distinct(probe$pid) == 2L, any(overlap))
piece_names <- tar_meta(
  names = partitions_patterns,
  fields = children,
  store = store
)$children[[1L]]
stopifnot(length(piece_names) == 3L)
piece_path <- file.path(store, "objects", piece_names[1L])
stopifnot(startsWith(
  normalizePath(piece_path, winslash = "/"),
  paste0(normalizePath(store, winslash = "/"), "/")
))
stopifnot(file.remove(piece_path))
outdated <- tar_outdated(script = script, store = store, reporter = "silent")
stopifnot(any(startsWith(outdated, "partitions_patterns")))
stopifnot(identical(
  list.files(file.path(store, "user"), recursive = TRUE),
  legacy_files
))
cat(
  "PASS: legacy result upgrade, native partitions, parallel IGS/ISS/ISPS/Mantel, partial Mantel wrapper, mapped outcomes, raw EEG regions, overlapping workers and missing-object invalidation.\n"
)
