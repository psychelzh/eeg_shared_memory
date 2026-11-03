# Optional real-cache check: three CCA components, one actual time window,
# and the production 9,999 inner Mantel permutations. No full pipeline run.
library(targets)
library(tidyverse)
tar_source()
group <- tar_read_raw("patterns_group_dynamic")
group <- filter(group, time_id == min(time_id))
yd <- tar_read_raw("pattern_shapes_alexnet_conv2d_4_9")
stopifnot(nrow(group) == 3L)
directory <- tempfile("real-mantel-partitions-")
bundle <- write_table_partitions(group, "cca_id", directory)
jobs <- make_partition_jobs(list(group = bundle))
sources <- c(group = "unused")
expr <- quote(extract_stats_mantel(calc_mantel(group, permute_dist(yd))))
fun <- eval(partition_callback(expr, sources))
plan <- make_permutation_plan(1L, 1L, seed = 71L)
seed <- plan$seed[plan$tar_rep == 1L]
reference <- with_permutation_seed(seed, fun(list(group = group), seed)) |>
  arrange(cca_id, time_id)
actual <- map(seq_len(nrow(jobs)), function(i) {
  run_partition_batch(jobs[i, ], plan, 1L, fun)
}) |>
  list_rbind() |>
  select(all_of(names(reference))) |>
  arrange(cca_id, time_id)
stopifnot(isTRUE(all.equal(reference, actual, tolerance = 1e-12)))
models <- calc_mantel(group, yd, seed = plan$seed[1L])
stopifnot(all(map_int(models$mantel, ~ length(.x$perm)) == 9999L))
cat(
  "PASS: real 150-stimulus GWS distances, three components and 9,999 inner Mantel permutations match the full-table reference.\n"
)
