library(tidyverse)
source("R/utils.R")
source("R/patterns.R")
source("R/behav.R")
source("R/prediction.R")
source("R/clusters.R")
source("R/mantel.R")
source("R/partitioning.R")
source("R/tar_partitions.R")

set.seed(47)
fixture <- crossing(
  cca_id = c(1L, 3L),
  subj_id = 1:8,
  time_id = seq(51L, 76L, 5L)
) |>
  mutate(pattern = map(subj_id, ~ dist(matrix(rnorm(144), 12))))
fixture$pattern[[2]][5] <- NA_real_
group <- fixture |>
  mutate(pattern = map(pattern, ~ dist(matrix(rnorm(144), 12))))
target_y <- dist(matrix(rnorm(144), 12))
target_z <- dist(matrix(rnorm(144), 12))
memory <- tibble(subj_id = 1:8, dprime = rnorm(8))
numeric_data <- corr_patterns_1(fixture, target_y, name = "iss")
directory <- tempfile("partition-check-")
bundles <- list(
  fixture = write_table_partitions(
    fixture,
    "cca_id",
    file.path(directory, "individual")
  ),
  group = write_table_partitions(
    group,
    "cca_id",
    file.path(directory, "group")
  ),
  numeric_data = write_table_partitions(
    numeric_data,
    "cca_id",
    file.path(directory, "numeric")
  )
)
jobs <- make_partition_jobs(bundles)
set.seed(71L)
expected_rng <- runif(n_distinct(fixture$time_id))
set.seed(71L)
native_rng <- run_partition_job(
  jobs[1L, ],
  function(.inputs, .replicate_seed) {
    tibble(value = runif(1))
  },
  seed = NULL
)
stopifnot(identical(native_rng$value, expected_rng))
plan <- make_permutation_plan(6L, 2L, 123L)
stopifnot(
  nrow(jobs) == 2L,
  !anyDuplicated(plan$seed),
  all(vapply(jobs$.files, length, integer(1)) == 3L)
)
row_seeds <- partition_row_seeds(fixture, 71L)
stopifnot(
  !anyDuplicated(row_seeds),
  identical(
    row_seeds[rev(seq_len(nrow(fixture)))],
    partition_row_seeds(fixture[rev(seq_len(nrow(fixture))), ], 71L)
  )
)
identical_inputs <- rep(list(target_y), nrow(fixture))
independent <- permute_pattern_rows(identical_inputs, fixture, 71L)
stopifnot(length(unique(map_chr(independent, digest::digest))) > 1L)
sources <- setNames(rep("unused", 3L), names(bundles))
inputs <- list(fixture = fixture, group = group, numeric_data = numeric_data)
order_result <- function(x) {
  x |>
    arrange(across(any_of(c(
      "cca_id",
      "region_id",
      "tar_rep",
      "time_id",
      "subj_id"
    ))))
}
compare <- function(expr, label, selected_jobs = jobs, values = inputs) {
  fun <- eval(partition_callback(expr, sources))
  expected <- map(seq_len(nrow(plan) - 1L) + 1L, function(i) {
    with_permutation_seed(plan$seed[i], fun(values, plan$seed[i])) |>
      mutate(tar_batch = plan$tar_batch[i], tar_rep = plan$tar_rep[i])
  }) |>
    list_rbind() |>
    order_result()
  actual <- map(seq_len(nrow(selected_jobs)), function(i) {
    map(1:3, ~ run_partition_batch(selected_jobs[i, ], plan, .x, fun)) |>
      list_rbind()
  }) |>
    list_rbind() |>
    order_result()
  stopifnot(isTRUE(all.equal(expected, actual, tolerance = 1e-12)))
  cat(
    "PASS:",
    label,
    "full-table vs streamed partitions under identical random inputs.\n"
  )
  actual
}
null_iss <- compare(
  quote(calc_stats_t(
    corr_patterns_1(fixture, permute_dist(target_y), name = "iss"),
    iss,
    alternative = "greater"
  )),
  "ISS/IWS shared target shuffle"
)
shared_reference <- with_permutation_seed(plan$seed[2L], {
  calc_stats_t(
    corr_patterns_1(fixture, permute_dist(target_y), name = "iss"),
    iss,
    alternative = "greater"
  )
}) |>
  arrange(cca_id, time_id)
shared_actual <- filter(null_iss, tar_rep == 1L) |>
  select(all_of(names(shared_reference))) |>
  arrange(cca_id, time_id)
stopifnot(isTRUE(all.equal(shared_reference, shared_actual, tolerance = 1e-12)))
compare(
  quote(calc_stats_t(
    corr_patterns(
      fixture,
      group |> mutate(pattern = map(pattern, permute_dist)),
      name = "igs"
    ),
    igs,
    alternative = "greater"
  )),
  "IGS independent subject/window shuffles"
) |>
  invisible()
compare(
  quote(
    fixture |>
      mutate(pattern = map(pattern, permute_dist)) |>
      calc_isps() |>
      summarise_isps()
  ),
  "ISPS with missing distances"
) |>
  invisible()
compare(
  quote(corr_mem(
    numeric_data,
    mutate(memory, subj_id = sample(subj_id)),
    col_pred = iss,
    alternative = "greater"
  )),
  "memory labels shared across components/time"
) |>
  invisible()
compare(
  quote(extract_stats_mantel(calc_mantel(
    group |> filter(subj_id == 1L),
    permute_dist(target_y),
    permutations = 19L
  ))),
  "Mantel shared outer target and independent inner row seeds"
) |>
  invisible()
compare(
  quote(extract_stats_mantel(calc_mantel_partial(
    group |> filter(subj_id == 1L),
    permute_dist(target_y),
    target_z,
    permutations = 19L
  ))),
  "partial Mantel"
) |>
  invisible()
compare(
  quote(
    regress_patterns_1(fixture, target_y) |> calc_isps() |> summarise_isps()
  ),
  "semantic control ISPS"
) |>
  invisible()
compare(
  quote(regress_patterns(fixture, group) |> calc_isps() |> summarise_isps()),
  "group-pattern control ISPS"
) |>
  invisible()
compare(
  quote(corr_patterns(
    regress_patterns_1(fixture, target_y),
    regress_patterns_1(group, target_y),
    name = "igs"
  )),
  "partial IGS"
) |>
  invisible()
compare(
  quote(corr_patterns(
    regress_patterns(fixture, group),
    regress_patterns_1r(target_y, group),
    name = "iss"
  )),
  "partial ISS"
) |>
  invisible()

# Full-series cluster inference is performed only after window results combine.
observed <- calc_stats_t(
  corr_patterns_1(fixture, target_y, name = "iss"),
  iss
) |>
  mutate(p.value = convert_p2_p1(p.value, statistic))
fun <- eval(partition_callback(
  quote(calc_stats_t(
    corr_patterns_1(
      fixture,
      permute_dist(target_y),
      name = "iss"
    ),
    iss,
    alternative = "greater"
  )),
  sources
))
rebatched <- make_permutation_plan(6L, 3L, 123L)
null_rebatched <- map(seq_len(nrow(jobs)), function(i) {
  map(1:2, ~ run_partition_batch(jobs[i, ], rebatched, .x, fun)) |> list_rbind()
}) |>
  list_rbind() |>
  order_result()
stopifnot(isTRUE(all.equal(
  select(null_iss, -tar_batch),
  select(null_rebatched, -tar_batch)
)))
stopifnot(isTRUE(all.equal(
  calc_clusters_stats(observed, null_iss),
  calc_clusters_stats(observed, null_rebatched)
)))

# Whole-period references retain the same shuffle at every dynamic time point.
whole_group <- filter(group, time_id == 51L) |> select(-time_id)
whole_bundles <- bundles
whole_bundles$group <- write_table_partitions(
  whole_group,
  "cca_id",
  file.path(directory, "whole")
)
whole_inputs <- inputs
whole_inputs$group <- whole_group
compare(
  quote(corr_patterns(
    fixture,
    group |> mutate(pattern = map(pattern, permute_dist)),
    name = "igs"
  )),
  "RT whole-period LOO reference",
  make_partition_jobs(whole_bundles),
  whole_inputs
) |>
  invisible()
cat(
  "PASS: rebatching, complete time clusters, shared and independent randomisation policies.\n"
)
