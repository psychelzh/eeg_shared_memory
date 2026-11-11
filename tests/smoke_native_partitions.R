# Opt-in native-cache alternative: small generated data and temporary projects.
library(targets)
library(tidyverse)
tar_source()

test_root <- tempfile("native-partitions-")
project_a <- file.path(test_root, "machine-a")
project_b <- file.path(test_root, "machine-b")
dir.create(project_a, recursive = TRUE)
fs::dir_copy("R", file.path(project_a, "R"))
stopifnot(file.copy(
  "tests/fixtures/native_partition_pipeline.R",
  file.path(project_a, "_targets.R")
))
sort_result <- function(x) arrange(x, cca_id, tar_rep, time_id)
withr::with_dir(project_a, {
  tar_make(reporter = "silent")
  partition_objects <- tar_read_raw("patterns")
  grouped_objects <- tar_read_raw("group")
  native_stats <- tar_read_raw("stats")
  plan <- tar_read_raw("plan")
  yd <- tar_read_raw("yd")
  probe <- map(native_stats, "probe") |> list_rbind()
  stopifnot(
    length(partition_objects) == 3L,
    all(map_lgl(partition_objects, ~ n_distinct(.x$cca_id) == 1L)),
    nrow(probe) == 12L,
    all(probe$rows_loaded == 32L),
    all(probe$pid != Sys.getpid()),
    length(list.files("_targets/user", recursive = TRUE)) == 0L
  )
  # The materialized native branches really have independent stored objects.
  meta <- tar_meta(fields = c(name, type, parent))
  pieces <- filter(meta, type == "branch", parent == "patterns")$name
  stopifnot(
    length(pieces) == 3L,
    all(file.exists(file.path("_targets/objects", pieces)))
  )
  inputs <- list(
    patterns = list_rbind(partition_objects),
    group = list_rbind(grouped_objects)
  )
  # Compare native execution with the current file-backed implementation and seeds.
  bundles <- imap(
    inputs,
    ~ write_table_partitions(
      .x,
      "cca_id",
      file.path("reference-bundles", .y)
    )
  )
  jobs <- make_partition_jobs(bundles)
  for (kind in c("iss", "igs")) {
    expr <- if (kind == "iss") {
      quote(calc_stats_t(
        corr_patterns_1(patterns, permute_dist(yd), name = "iss"),
        iss,
        alternative = "greater"
      ))
    } else {
      quote(calc_stats_t(
        corr_patterns(
          patterns,
          mutate(group, pattern = map(pattern, permute_dist)),
          name = "igs"
        ),
        igs,
        alternative = "greater"
      ))
    }
    fun <- eval(partition_callback(
      expr,
      c(patterns = "unused", group = "unused")
    ))
    reference <- crossing(job = seq_len(nrow(jobs)), batch = 1:4) |>
      pmap(function(job, batch) {
        run_partition_batch(jobs[job, ], plan, batch, fun)
      }) |>
      list_rbind() |>
      sort_result()
    actual <- map(native_stats, kind) |> list_rbind() |> sort_result()
    stopifnot(isTRUE(all.equal(reference, actual, tolerance = 1e-12)))
  }
})

fs::dir_copy(project_a, project_b)
withr::with_dir(project_b, {
  stopifnot(length(tar_outdated(reporter = "silent")) == 0L)
  tar_make(reporter = "silent")
  stopifnot(identical(native_stats, tar_read_raw("stats")))
  # Standard tar_prune owns this cache: obsolete native objects need no directory scan.
  script <- readLines("_targets.R")
  script <- script[!grepl("tar_target\\(obsolete_probe", script)]
  # Remove the trailing comma belonging to stats after dropping the last target.
  script[length(script) - 1L] <- sub("\\),$", ")", script[length(script) - 1L])
  writeLines(script, "_targets.R")
  stopifnot(file.exists("_targets/objects/obsolete_probe"))
  tar_prune(cloud = FALSE)
  stopifnot(
    !file.exists("_targets/objects/obsolete_probe"),
    !"obsolete_probe" %in% tar_meta(fields = name)$name,
    identical(native_stats, tar_read_raw("stats"))
  )
})
cat(
  "PASS: native per-CCA storage, worker subset reads, crossed batches, ISS/IGS numerical equivalence, copied-store cache reuse and native pruning.\n"
)
