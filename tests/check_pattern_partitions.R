library(tidyverse)
source("R/utils.R")
source("R/patterns.R")
source("R/io.R")
source("R/partitioning.R")
source("R/tar_partitions.R")
set.seed(61)
signals <- crossing(
  cca_id = 1:2,
  subj_id = 1:4,
  time_id = 1:61,
  trial_id = 1:8
) |>
  mutate(y = rnorm(n()))
sort_result <- function(x) {
  x |> arrange(across(any_of(c("region_id", "subj_id", "cca_id", "time_id"))))
}
for (fun in list(
  calc_indiv_pattern,
  calc_group_pattern,
  calc_indiv_pattern_dynamic,
  calc_group_pattern_dynamic
)) {
  reference <- fun(signals) |> sort_result()
  partitioned <- split(signals, signals$cca_id) |>
    map(fun) |>
    list_rbind() |>
    sort_result()
  stopifnot(isTRUE(all.equal(reference, partitioned, tolerance = 1e-12)))
}
directory <- tempfile("region-partitions-")
dir.create(directory)
data <- array(rnorm(12 * 61 * 8 * 6), c(12, 61, 8, 6))
file <- file.path(directory, "raw.qs2")
qs2::qs_save(data, file, nthreads = 1L)
regions <- setNames(
  lapply(1:6, function(i) as.integer(seq_len(12) %in% c(2 * i - 1, 2 * i))),
  paste0("region", 1:6)
)
reference <- read_eeg_regions(file, regions)
bundle <- write_eeg_region_partitions(
  file,
  regions,
  file.path(directory, "pieces")
)
jobs <- make_partition_jobs(list(eeg = bundle))
for (fun in list(
  calc_indiv_pattern_region,
  calc_indiv_pattern_dynamic_region
)) {
  expected <- reference |>
    mutate(data = map(data, fun)) |>
    unnest(data) |>
    sort_result()
  actual <- map(seq_len(nrow(jobs)), function(i) {
    run_partition_job(
      jobs[i, ],
      function(.inputs, .replicate_seed) fun(.inputs$eeg),
      seed = 13L
    )
  }) |>
    list_rbind() |>
    sort_result()
  stopifnot(isTRUE(all.equal(expected, actual, tolerance = 1e-12)))
}
# Regions with one channel remain four-dimensional when written.
regions$region1 <- as.integer(seq_len(12) == 1L)
single <- write_eeg_region_partitions(
  file,
  regions,
  file.path(directory, "single")
)
stopifnot(identical(
  dim(qs2::qs_read(single[2], nthreads = 1L)),
  c(1L, 61L, 8L, 6L)
))
cat(
  "PASS: CCA pattern construction, raw-region dynamic/whole calculations, and array dimensions.\n"
)
