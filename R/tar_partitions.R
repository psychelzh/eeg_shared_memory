partition_sources_for <- function(expr, sources = NULL) {
  if (!is.null(sources)) {
    return(sources)
  }
  registry <- get0(
    "partition_sources",
    envir = parent.frame(),
    inherits = TRUE,
    ifnotfound = character()
  )
  registry[intersect(names(registry), all.names(expr))]
}

# Only known randomisation operations are rewritten. Shared target/behaviour
# shuffles still use the replicate seed; individual pattern rows use key seeds.
partition_expression <- function(expr, sources) {
  if (is.symbol(expr)) {
    if (as.character(expr) %in% names(sources)) {
      return(rlang::call2("$", quote(.inputs), expr))
    }
    return(expr)
  }
  if (!is.call(expr)) {
    return(expr)
  }
  if (identical(expr[[1L]], quote(`$`)) || identical(expr[[1L]], quote(`@`))) {
    expr[[2L]] <- partition_expression(expr[[2L]], sources)
    return(expr)
  }
  if (
    (identical(expr[[1L]], quote(map)) ||
      identical(expr[[1L]], quote(purrr::map))) &&
      length(expr) >= 3L &&
      identical(expr[[2L]], quote(pattern)) &&
      "permute_dist" %in% all.names(expr[[3L]])
  ) {
    if (length(expr) != 3L || !identical(expr[[3L]], quote(permute_dist))) {
      stop(
        "Use map(pattern, permute_dist) for independent keyed pattern permutations."
      )
    }
    return(quote(permute_pattern_rows(
      pattern,
      dplyr::pick(dplyr::everything()),
      .replicate_seed
    )))
  }
  expr[] <- lapply(as.list(expr), partition_expression, sources = sources)
  if (
    is.symbol(expr[[1L]]) &&
      as.character(expr[[1L]]) %in% c("calc_mantel", "calc_mantel_partial")
  ) {
    expr$seed <- quote(.replicate_seed)
  }
  expr
}

partition_callback <- function(expr, sources) {
  as.call(list(
    quote(`function`),
    as.pairlist(alist(.inputs = , .replicate_seed = )),
    partition_expression(expr, sources)
  ))
}

partition_jobs_target <- function(name, sources) {
  catalogs <- rlang::call2(
    "list",
    !!!setNames(
      lapply(unname(sources), function(source) {
        as.name(paste0(source, "_keys"))
      }),
      names(sources)
    )
  )
  targets::tar_target_raw(
    name,
    rlang::call2("make_native_partition_jobs", catalogs),
    iteration = "vector",
    retrieval = "worker",
    memory = "transient"
  )
}

partition_inputs_expression <- function(jobs, sources) {
  rlang::call2(
    "native_partition_inputs",
    as.name(jobs),
    rlang::call2(
      "list",
      !!!setNames(lapply(unname(sources), as.name), names(sources))
    )
  )
}

partition_map_pattern <- function(jobs, sources) {
  rlang::call2("map", !!!lapply(c(jobs, unique(unname(sources))), as.name))
}

tar_native_partitions <- function(
  name,
  data,
  by = "cca_id",
  resources = targets::tar_option_get("resources")
) {
  keys <- paste0(name, "_keys")
  list(
    targets::tar_target_raw(
      keys,
      bquote(native_partition_keys(.(substitute(data)), .(by))),
      iteration = "vector",
      resources = resources,
      retrieval = "worker",
      memory = "transient",
      garbage_collection = TRUE
    ),
    targets::tar_target_raw(
      name,
      bquote(native_partition_piece(.(substitute(data)), .(as.name(keys)))),
      pattern = rlang::call2("map", as.name(keys)),
      iteration = "list",
      resources = resources,
      retrieval = "worker",
      memory = "transient",
      garbage_collection = TRUE
    )
  )
}

tar_native_eeg_partitions <- function(name, file, channel_regions, resources) {
  keys <- paste0(name, "_keys")
  list(
    targets::tar_target_raw(
      keys,
      quote(tibble::tibble(region_id = seq_len(6L))),
      iteration = "vector",
      resources = resources
    ),
    targets::tar_target_raw(
      name,
      bquote(native_eeg_partition_piece(
        .(substitute(file)),
        .(substitute(channel_regions)),
        .(as.name(keys))
      )),
      pattern = rlang::call2("map", as.name(keys)),
      iteration = "list",
      resources = resources,
      retrieval = "worker",
      memory = "transient",
      garbage_collection = TRUE
    )
  )
}

tar_partition_registry <- function(
  sources,
  by = list(),
  default_by = "cca_id",
  resources = targets::tar_option_get("resources")
) {
  lapply(names(sources), function(source) {
    keys <- by[[source]]
    if (is.null(keys)) {
      keys <- default_by
    }
    eval(bquote(tar_native_partitions(
      .(unname(sources[source])),
      .(as.name(source)),
      by = keys,
      resources = resources
    )))
  })
}

tar_partition_target <- function(name, command, sources = NULL, stream = TRUE) {
  expr <- substitute(command)
  sources <- partition_sources_for(expr, sources)
  jobs <- paste0("jobs_", name)
  seed <- paste0("seed_", name)
  stochastic <- any(
    all.names(expr) %in% c("calc_mantel", "calc_mantel_partial")
  )
  seed_expr <- if (stochastic) as.name(seed) else NULL
  list(
    partition_jobs_target(jobs, sources),
    if (stochastic) {
      targets::tar_target_raw(
        seed,
        quote(make_permutation_plan(1L, 1L, targets::tar_seed_get())$seed[1L])
      )
    },
    targets::tar_target_raw(
      name,
      bquote(run_partition_job(
        .(as.name(jobs)),
        .(partition_callback(expr, sources)),
        .(seed_expr),
        stream = .(stream),
        inputs = .(partition_inputs_expression(jobs, sources))
      )),
      pattern = partition_map_pattern(jobs, sources),
      memory = "transient",
      retrieval = "worker",
      garbage_collection = TRUE
    )
  )
}

partition_batch_targets <- function(name, expr, sources, reps, batch_size) {
  jobs <- paste0("jobs_", name)
  plan <- paste0("plan_", name)
  batches <- paste0("batches_", name)
  list(
    partition_jobs_target(jobs, sources),
    targets::tar_target_raw(
      plan,
      bquote(make_permutation_plan(
        .(reps),
        .(batch_size),
        targets::tar_seed_get()
      ))
    ),
    targets::tar_target_raw(
      batches,
      bquote(sort(unique(.(as.name(plan))$tar_batch[-1L])))
    ),
    targets::tar_target_raw(
      name,
      bquote(run_partition_batch(
        .(as.name(jobs)),
        .(as.name(plan)),
        .(as.name(batches)),
        .(partition_callback(expr, sources)),
        inputs = .(partition_inputs_expression(jobs, sources))
      )),
      pattern = rlang::call2(
        "cross",
        partition_map_pattern(jobs, sources),
        as.name(batches)
      ),
      memory = "transient",
      retrieval = "worker",
      garbage_collection = TRUE
    )
  )
}

tar_partition_replicates <- function(
  name,
  command,
  reps = 1000L,
  batch_size = 10L,
  sources = NULL
) {
  expr <- substitute(command)
  sources <- partition_sources_for(expr, sources)
  partition_batch_targets(name, expr, sources, reps, batch_size)
}

tar_cluster_permutation <- function(
  name,
  stats_expr,
  stats_perm_expr,
  data_expr = NULL,
  data_perm_expr = NULL,
  clusters_stats_expr = NULL,
  stats_name = NULL,
  data_name = NULL,
  clusters_stats_name = NULL,
  reps = 10L,
  batches = 100L,
  sources = NULL,
  batch_size = NULL
) {
  if (!missing(name)) {
    stats_name <- paste0("stats_", name)
    data_name <- paste0("data_", name)
    clusters_stats_name <- paste0("clusters_stats_", name)
  }
  observed <- substitute(stats_expr)
  permuted <- substitute(stats_perm_expr)
  data <- substitute(data_expr)
  data_permuted <- substitute(data_perm_expr)
  sources <- partition_sources_for(
    as.call(list(quote(list), data, data_permuted, observed, permuted)),
    sources
  )
  stats_permuted_name <- paste0(stats_name, "_permuted")
  if (!is.null(data)) {
    observed <- targets::tar_tidy_eval(
      as.expression(observed),
      envir = list(.x = as.name(data_name)),
      tidy_eval = TRUE
    )
    permuted <- targets::tar_tidy_eval(
      as.expression(permuted),
      envir = list(.x = quote(.permuted)),
      tidy_eval = TRUE
    )
    permuted <- bquote(local({
      .permuted <- .(data_permuted)
      .(permuted)
    }))
  }
  if (is.null(batch_size)) {
    batch_size <- if (
      any(all.names(permuted) %in% c("calc_mantel", "calc_mantel_partial"))
    ) {
      min(reps, 2L)
    } else {
      reps
    }
  }
  targets_perm <- partition_batch_targets(
    stats_permuted_name,
    permuted,
    sources,
    reps * batches,
    batch_size
  )
  clusters <- substitute(clusters_stats_expr)
  if (is.null(clusters)) {
    clusters <- rlang::call2(
      "calc_clusters_stats",
      as.name(stats_name),
      as.name(stats_permuted_name)
    )
  } else {
    clusters <- targets::tar_tidy_eval(
      as.expression(clusters),
      envir = list(.x = as.name(stats_name), .y = as.name(stats_permuted_name)),
      tidy_eval = TRUE
    )
  }
  observed_seed <- if (
    any(all.names(data) %in% c("calc_mantel", "calc_mantel_partial"))
  ) {
    bquote(.(as.name(paste0("plan_", stats_permuted_name)))$seed[1L])
  } else {
    NULL
  }
  list(
    targets_perm,
    if (!is.null(data)) {
      targets::tar_target_raw(
        data_name,
        bquote(run_partition_job(
          .(as.name(paste0("jobs_", stats_permuted_name))),
          .(partition_callback(data, sources)),
          .(observed_seed),
          inputs = .(partition_inputs_expression(
            paste0("jobs_", stats_permuted_name),
            sources
          ))
        )),
        pattern = partition_map_pattern(
          paste0("jobs_", stats_permuted_name),
          sources
        ),
        memory = "transient",
        retrieval = "worker",
        garbage_collection = TRUE
      )
    },
    targets::tar_target_raw(stats_name, observed, retrieval = "worker"),
    targets::tar_target_raw(clusters_stats_name, clusters, retrieval = "worker")
  )
}
