# Isolated proof of native, persistent CCA partitions. No study inputs are read.
library(targets)
library(tidyverse)
tar_source()
tar_option_set(
  packages = "tidyverse",
  format = "qs",
  memory = "transient",
  retrieval = "worker",
  controller = crew::crew_controller_local(name = "native", workers = 2),
  resources = tar_resources(crew = tar_resources_crew(controller = "native"))
)

native_piece <- function(cca_id) {
  with_permutation_seed(100L + cca_id, {
    crossing(subj_id = 1:8, time_id = seq(51L, 66L, 5L)) |>
      mutate(
        cca_id = cca_id,
        pattern = map(subj_id, ~ dist(matrix(rnorm(144), 12)))
      )
  })
}

native_batch <- function(patterns, group, yd, plan, batch) {
  stopifnot(
    n_distinct(patterns$cca_id) == 1L,
    identical(unique(patterns$cca_id), unique(group$cca_id))
  )
  inputs <- list(patterns = patterns, group = group)
  sources <- c(patterns = "unused", group = "unused")
  expressions <- list(
    iss = quote(calc_stats_t(
      corr_patterns_1(patterns, permute_dist(yd), name = "iss"),
      iss,
      alternative = "greater"
    )),
    igs = quote(calc_stats_t(
      corr_patterns(
        patterns,
        mutate(group, pattern = map(pattern, permute_dist)),
        name = "igs"
      ),
      igs,
      alternative = "greater"
    ))
  )
  job <- tibble(.key = list(tibble(cca_id = unique(patterns$cca_id))))
  selected <- plan[plan$tar_batch == batch & plan$tar_rep > 0L, ]
  result <- lapply(expressions, function(expr) {
    fun <- eval(partition_callback(expr, sources))
    map(seq_len(nrow(selected)), function(i) {
      run_partition_job(job, fun, selected$seed[i], inputs = inputs) |>
        mutate(tar_batch = batch, tar_rep = selected$tar_rep[i])
    }) |>
      list_rbind()
  })
  result$probe <- tibble(
    cca_id = unique(patterns$cca_id),
    batch = batch,
    rows_loaded = nrow(patterns),
    pid = Sys.getpid()
  )
  result
}

list(
  tar_target(cca_id, 1:3),
  tar_target(
    patterns,
    native_piece(cca_id),
    pattern = map(cca_id),
    iteration = "list"
  ),
  tar_target(
    group,
    with_permutation_seed(200L + unique(patterns$cca_id), {
      mutate(patterns, pattern = map(pattern, ~ dist(matrix(rnorm(144), 12))))
    }),
    pattern = map(patterns),
    iteration = "list"
  ),
  tar_target(yd, with_permutation_seed(71L, dist(matrix(rnorm(144), 12)))),
  tar_target(plan, make_permutation_plan(8L, 2L, 72L)),
  tar_target(batch, 1:4),
  tar_target(
    stats,
    native_batch(patterns, group, yd, plan, batch),
    pattern = cross(map(patterns, group), batch),
    iteration = "list"
  ),
  tar_target(obsolete_probe, 17L)
)
