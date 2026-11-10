# Exercise relocation in two isolated projects, never the completed study cache.
library(tidyverse)
library(targets)
source("R/partitioning.R")

test_root <- tempfile("partition-portability-")
dir.create(test_root)
project_a <- file.path(test_root, "machine-a", "study")
project_b <- file.path(test_root, "machine-b", "study")
dir.create(file.path(project_a, "R"), recursive = TRUE)
stopifnot(file.copy(
  "R/partitioning.R",
  file.path(project_a, "R", "partitioning.R")
))
fixture_data <- tibble(cca_id = rep(1:2, each = 3), value = seq_len(6))
saveRDS(fixture_data, file.path(project_a, "input.rds"), version = 2L)
writeLines(
  c(
    "library(targets)",
    "library(tidyverse)",
    "source('R/partitioning.R', local = TRUE)",
    "tar_option_set(packages = 'tidyverse', controller = NULL)",
    "list(",
    "  tar_target(input_file, 'input.rds', format = 'file'),",
    "  tar_target(input, readRDS(input_file)),",
    "  tar_target(bundle, write_table_partitions(input, 'cca_id',",
    "    '_targets/user/partitions/input'), format = 'file'),",
    "  tar_target(jobs, make_partition_jobs(list(data = bundle))),",
    "  tar_target(result, {",
    "    write('ran', 'executions.txt', append = TRUE)",
    "    run_partition_job(jobs, function(.inputs, .replicate_seed) {",
    "      tibble(total = sum(.inputs$data$value))",
    "    }, seed = 71L)",
    "  }, pattern = map(jobs))",
    ")"
  ),
  file.path(project_a, "_targets.R")
)

withr::with_dir(project_a, {
  tar_make(callr_function = NULL, reporter = "silent", envir = new.env())
  original <- tar_read_raw("result")
  bundle_paths <- tar_read_raw("bundle")
  cached_jobs <- tar_read_raw("jobs")
  stopifnot(
    !any(grepl("^([A-Za-z]:|/)", bundle_paths)),
    !any(grepl("^([A-Za-z]:|/)", unlist(cached_jobs$.files))),
    identical(readRDS(bundle_paths[1L])$.file, c("part-1.qs2", "part-2.qs2"))
  )
  # Identical preparation must preserve paths, bytes and existing modification times.
  before <- file.info(bundle_paths)$mtime
  repeated <- write_table_partitions(
    fixture_data,
    "cca_id",
    "_targets/user/partitions/input"
  )
  stopifnot(
    identical(bundle_paths, repeated),
    identical(before, file.info(bundle_paths)$mtime)
  )
  original_branches <- tar_meta(names = result, fields = children)$children[[
    1L
  ]]
  # Missing pieces are restored within the same content-addressed generation.
  stopifnot(file.remove(bundle_paths[2L]))
  repaired <- write_table_partitions(
    fixture_data,
    "cca_id",
    "_targets/user/partitions/input"
  )
  stopifnot(identical(bundle_paths, repaired), all(file.exists(bundle_paths)))
})

dir.create(dirname(project_b), recursive = TRUE)
fs::dir_copy(project_a, project_b)
withr::with_dir(project_b, {
  # The old root still exists, so accidental reads from it cannot make this test pass.
  stopifnot(
    length(tar_outdated(callr_function = NULL, reporter = "silent")) == 0L
  )
  before <- readLines("executions.txt")
  tar_make(callr_function = NULL, reporter = "silent", envir = new.env())
  stopifnot(
    identical(before, readLines("executions.txt")),
    identical(original, tar_read_raw("result")),
    identical(
      original_branches,
      tar_meta(names = result, fields = children)$children[[1L]]
    )
  )
  # Legacy manifest entries may contain either Windows or Linux absolute paths.
  manifest <- readRDS(bundle_paths[1L])
  for (prefix in c("Z:/retired-device/study/", "/retired-server/study/")) {
    legacy <- manifest
    legacy$.file <- paste0(prefix, dirname(bundle_paths[1L]), "/", legacy$.file)
    legacy_manifest <- file.path(
      dirname(bundle_paths[1L]),
      "legacy-manifest.rds"
    )
    saveRDS(legacy, legacy_manifest, version = 2L)
    migrated_jobs <- make_partition_jobs(list(
      data = c(legacy_manifest, bundle_paths[-1L])
    ))
    stopifnot(identical(cached_jobs, migrated_jobs))
    legacy_job <- cached_jobs[1L, ]
    legacy_job$.files[[1L]] <- setNames(
      paste0(prefix, legacy_job$.files[[1L]]),
      names(legacy_job$.files[[1L]])
    )
    stopifnot(identical(
      read_partition_inputs(cached_jobs[1L, ]),
      read_partition_inputs(legacy_job)
    ))
  }
  # An actual input change must change the bundle identity and recompute results.
  changed <- fixture_data
  changed$value[1L] <- 100L
  saveRDS(changed, "input.rds", version = 2L)
  tar_make(callr_function = NULL, reporter = "silent", envir = new.env())
  stopifnot(
    !identical(bundle_paths, tar_read_raw("bundle")),
    !identical(original, tar_read_raw("result")),
    length(readLines("executions.txt")) > length(before),
    sum(tar_read_raw("result")$total) == sum(changed$value)
  )
})
cat(
  "PASS: portable file/job paths, content-stable bundles, repair, legacy path reads, copied-store cache hits and genuine-change invalidation.\n"
)
