# Production factories, alignment failures and relocation in temporary stores.
library(targets)
library(tidyverse)
source("R/partitioning.R")
source("R/tar_partitions.R")

fixture_data <- tibble(cca_id = rep(1:2, each = 3), value = seq_len(6))
catalog <- native_partition_keys(fixture_data[6:1, ], "cca_id")
stopifnot(identical(catalog$cca_id, 1:2))
key_jobs <- make_native_partition_jobs(list(a = catalog, b = catalog))
stopifnot(
  nrow(key_jobs) == 2L,
  inherits(
    try(
      make_native_partition_jobs(list(a = catalog, b = catalog[2:1, ])),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(
      native_partition_inputs(
        key_jobs[1L, ],
        list(a = native_partition_piece(fixture_data, catalog[2L, ]))
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  identical(make_native_partition_jobs(list())$.key, list(tibble()))
)

test_root <- tempfile("native-factory-check-")
project_a <- file.path(test_root, "machine-a")
project_b <- file.path(test_root, "machine-b")
dir.create(file.path(project_a, "R"), recursive = TRUE)
stopifnot(all(file.copy(
  c("R/partitioning.R", "R/tar_partitions.R"),
  file.path(project_a, "R")
)))
saveRDS(fixture_data, file.path(project_a, "input.rds"))
writeLines(
  c(
    "library(targets)",
    "library(tidyverse)",
    "source('R/partitioning.R', local = TRUE)",
    "source('R/tar_partitions.R', local = TRUE)",
    "tar_option_set(packages = 'tidyverse', format = 'qs', controller = NULL)",
    "partition_sources <- c(input = 'partitions_input')",
    "list(",
    "  tar_target(input_file, 'input.rds', format = 'file'),",
    "  tar_target(input, readRDS(input_file)),",
    "  tar_partition_registry(partition_sources),",
    "  tar_partition_target('result', {",
    "    write('ran', 'executions.txt', append = TRUE)",
    "    tibble(total = sum(input$value))",
    "  })",
    ")"
  ),
  file.path(project_a, "_targets.R")
)

withr::with_dir(project_a, {
  tar_make(callr_function = NULL, reporter = "silent", envir = new.env())
  cached_result <- tar_read_raw("result")
  piece_names <- tar_meta(
    names = partitions_input,
    fields = children
  )$children[[1L]]
  pieces <- tar_read_raw("partitions_input")
  cached_jobs <- tar_read_raw("jobs_result")
  stopifnot(
    length(pieces) == 2L,
    all(vapply(pieces, function(x) nrow(x$data), integer(1)) == 3L),
    identical(names(cached_jobs), ".key"),
    !grepl(
      project_a,
      paste(capture.output(dput(pieces)), collapse = ""),
      fixed = TRUE
    ),
    length(list.files("_targets/user", recursive = TRUE)) == 0L
  )
  before <- readLines("executions.txt")
  piece_path <- file.path("_targets", "objects", piece_names[1L])
  stopifnot(file.remove(piece_path))
  stopifnot(any(startsWith(
    tar_outdated(callr_function = NULL, reporter = "silent"),
    "partitions_input"
  )))
  tar_make(callr_function = NULL, reporter = "silent", envir = new.env())
  stopifnot(
    file.exists(piece_path),
    identical(before, readLines("executions.txt"))
  )
})
fs::dir_copy(project_a, project_b)
withr::with_dir(project_b, {
  stopifnot(
    length(tar_outdated(callr_function = NULL, reporter = "silent")) == 0L
  )
  tar_make(callr_function = NULL, reporter = "silent", envir = new.env())
  stopifnot(
    identical(before, readLines("executions.txt")),
    identical(cached_result, tar_read_raw("result"))
  )
  changed <- fixture_data
  changed$value[1L] <- 100L
  saveRDS(changed, "input.rds")
  tar_make(callr_function = NULL, reporter = "silent", envir = new.env())
  stopifnot(
    sum(tar_read_raw("result")$total) == sum(changed$value),
    length(readLines("executions.txt")) > length(before)
  )
})
cat(
  "PASS: production native factories, ordered-key guards, no file payloads, missing-object repair, relocation cache hits and true-change invalidation.\n"
)
