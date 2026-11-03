# Explicit opt-in integration experiment: Rscript tests/smoke_worker_pools.R
# Uses its own temporary store, two workers, and tiny generated distances.
library(targets)
library(tidyverse)
source("R/mantel.R")
source("R/clusters.R")
source("R/patterns.R")
source("R/utils.R")

root <- normalizePath(".", winslash = "/")
script <- tempfile(fileext = ".R")
store <- tempfile(pattern = "eeg-pool-smoke-")
lines <- c(
  'library(targets)',
  sprintf('source(%s)', encodeString(paste0(root, "/R/tar_factories.R"), quote = '"')),
  sprintf('source(%s)', encodeString(paste0(root, "/R/mantel.R"), quote = '"')),
  sprintf('source(%s)', encodeString(paste0(root, "/R/clusters.R"), quote = '"')),
  'tar_option_set(packages = "tidyverse", format = "rds",',
  '  resources = tar_resources(crew = tar_resources_crew(controller = "light")),',
  '  controller = crew::crew_controller_group(',
  '    crew::crew_controller_local(name = "light", workers = 1),',
  '    crew::crew_controller_local(name = "heavy", workers = 1)))',
  'heavy <- tar_resources(crew = tar_resources_crew(controller = "heavy"))',
  'list(',
  '  tar_target(fixture, {set.seed(17); tibble(cca_id = 1L, time_id = 1:3,',
  '    pattern = replicate(3, dist(matrix(rnorm(60), 12)), simplify = FALSE))}),',
  '  tar_target(ydis, {set.seed(18); dist(matrix(rnorm(60), 12))}),',
  '  tar_cluster_permutation("smoke",',
  '    data_expr = local({set.seed(19); calc_mantel(fixture, ydis)}),',
  '    data_perm_expr = local({set.seed(29); calc_mantel(fixture, ydis)}),',
  '    stats_expr = extract_stats_mantel(!!.x),',
  '    stats_perm_expr = extract_stats_mantel(!!.x),',
  '    batches = 2, reps = 2, data_resources = heavy),',
  '  tar_target(light_probe, {start <- as.numeric(Sys.time()); Sys.sleep(4);',
  '    list(start = start, end = as.numeric(Sys.time()), pid = Sys.getpid())}),',
  '  tar_target(heavy_probe, {start <- as.numeric(Sys.time()); Sys.sleep(4);',
  '    list(start = start, end = as.numeric(Sys.time()), pid = Sys.getpid())}, resources = heavy)',
  ')'
)
writeLines(lines, script)
tar_make(script = script, store = store, reporter = "balanced")
fixture <- tar_read_raw("fixture", store = store)
ydis <- tar_read_raw("ydis", store = store)
set.seed(29)
full <- calc_mantel(fixture, ydis)
expected <- extract_stats_mantel(full)
permuted <- tar_read_raw("stats_smoke_permuted", store = store)
actual <- select(permuted, all_of(names(expected)))
stopifnot(nrow(actual) == 4L * nrow(expected),
          isTRUE(all.equal(actual, bind_rows(rep(list(expected), 4)))))
observed <- tar_read_raw("stats_smoke", store = store)
clusters_expected <- calc_clusters_stats(observed, permuted)
clusters_actual <- tar_read_raw("clusters_stats_smoke", store = store)
stopifnot(isTRUE(all.equal(clusters_actual, clusters_expected)))
light <- tar_read_raw("light_probe", store = store)
heavy <- tar_read_raw("heavy_probe", store = store)
stopifnot(light$pid != heavy$pid,
          max(light$start, heavy$start) < min(light$end, heavy$end))
cat(sprintf("PASS: 4 fused Mantel replicates and final clusters match reference; worker pools overlap.\nFull model %.1f KB -> statistics %.1f KB (%.1f%% smaller).\n",
            as.numeric(object.size(full)) / 1024,
            as.numeric(object.size(expected)) / 1024,
            100 * (1 - as.numeric(object.size(expected)) / as.numeric(object.size(full)))))

# Check the ISPS fusion across replicates, using the same random input in both
# paths. Replicate identifiers must survive summarisation unchanged.
set.seed(31)
individual <- crossing(subj_id = 1:20, cca_id = 1L, time_id = 1:3) |>
  mutate(pattern = map(subj_id, ~ dist(matrix(rnorm(60), 12))))
make_isps <- function() individual |>
  mutate(pattern = map(pattern, ~ .x[sample(length(.x))])) |>
  calc_isps()
set.seed(32)
full_isps <- map(1:4, ~ mutate(make_isps(), tar_rep = .x)) |> list_rbind()
reference_isps <- summarise_isps(full_isps)
set.seed(32)
fused_isps <- map(1:4, ~ mutate(summarise_isps(make_isps()), tar_rep = .x)) |>
  list_rbind()
stopifnot(isTRUE(all.equal(reference_isps, select(fused_isps, all_of(names(reference_isps))))))
cat(sprintf("PASS: ISPS summaries preserve all 4 replicates; %.1f KB -> %.1f KB.\n",
            as.numeric(object.size(full_isps)) / 1024,
            as.numeric(object.size(fused_isps)) / 1024))
