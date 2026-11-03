# A bundle tracks its manifest and every data file as one generated file target.
# Individual jobs carry only keys and file paths; workers load their own pieces.
write_partition_bundle <- function(keys, get_piece, directory) {
  stopifnot(nrow(keys) > 0L, !anyNA(keys), !anyDuplicated(keys))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  generation <- tempfile("generation-", tmpdir = directory)
  dir.create(generation)
  files <- vapply(
    seq_len(nrow(keys)),
    function(i) {
      file <- file.path(generation, paste0("part-", i, ".qs2"))
      qs2::qs_save(get_piece(i), file, nthreads = 1L)
      normalizePath(file, winslash = "/", mustWork = TRUE)
    },
    character(1)
  )
  manifest <- file.path(generation, "manifest.rds")
  saveRDS(dplyr::mutate(keys, .file = files), manifest)
  c(normalizePath(manifest, winslash = "/", mustWork = TRUE), files)
}

write_table_partitions <- function(data, by, directory) {
  stopifnot(is.data.frame(data), all(by %in% names(data)))
  groups <- vctrs::vec_group_loc(data[by])
  ordering <- do.call(order, groups$key)
  write_partition_bundle(
    groups$key[ordering, , drop = FALSE],
    function(i) {
      data[groups$loc[[ordering[i]]], , drop = FALSE]
    },
    directory
  )
}

write_eeg_region_partitions <- function(file, channel_regions, directory) {
  data <- qs2::qs_read(file, nthreads = 1L)
  keys <- tibble::tibble(region_id = seq_len(6L))
  write_partition_bundle(
    keys,
    function(i) {
      data[channel_regions[[paste0("region", i)]] != 0, , , , drop = FALSE]
    },
    directory
  )
}

make_partition_jobs <- function(bundles) {
  if (!length(bundles)) {
    return(tibble::tibble(
      .key = list(tibble::tibble()),
      .files = list(character())
    ))
  }
  manifests <- lapply(bundles, function(paths) readRDS(paths[1L]))
  by <- setdiff(names(manifests[[1L]]), ".file")
  stopifnot(all(vapply(
    manifests,
    function(x) {
      setequal(setdiff(names(x), ".file"), by) && !anyDuplicated(x[by])
    },
    logical(1)
  )))
  joined <- manifests[[1L]][by]
  for (name in names(manifests)) {
    part <- manifests[[name]]
    names(part)[names(part) == ".file"] <- name
    joined <- dplyr::full_join(
      joined,
      part,
      by = by,
      relationship = "one-to-one"
    )
  }
  stopifnot(!anyNA(joined))
  tibble::tibble(
    .key = lapply(seq_len(nrow(joined)), function(i) {
      joined[i, by, drop = FALSE]
    }),
    .files = lapply(seq_len(nrow(joined)), function(i) {
      unlist(joined[i, names(bundles), drop = FALSE], use.names = TRUE)
    })
  )
}

with_permutation_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  withr::with_seed(
    seed,
    code,
    .rng_kind = "Mersenne-Twister",
    .rng_sample_kind = "Rejection"
  )
}

make_permutation_plan <- function(reps, batch_size, seed) {
  stopifnot(
    reps > 0L,
    reps == as.integer(reps),
    batch_size > 0L,
    batch_size == as.integer(batch_size)
  )
  tibble::tibble(
    tar_rep = seq.int(0L, reps),
    tar_batch = c(0L, (seq_len(reps) - 1L) %/% batch_size + 1L),
    seed = with_permutation_seed(
      seed,
      sample.int(.Machine$integer.max, reps + 1L)
    )
  )
}

# Row keys keep independent subject/time randomisations invariant to partition
# boundaries, row order, batch size and which worker executes a task.
partition_row_seeds <- function(data, seed) {
  keys <- intersect(
    c("cca_id", "region_id", "subj_id", "time_id", "pair"),
    names(data)
  )
  if (!length(keys)) {
    stopifnot(nrow(data) == 1L)
    return(as.integer(seed))
  }
  labels <- lapply(keys, function(key) {
    paste0(key, "=", encodeString(as.character(data[[key]]), quote = '"'))
  })
  label <- do.call(paste, c(labels, sep = ";"))
  stopifnot(!anyDuplicated(label))
  digest::digest2int(label, seed = as.integer(seed))
}

permute_pattern_rows <- function(pattern, data, seed) {
  purrr::map2(pattern, partition_row_seeds(data, seed), function(x, row_seed) {
    with_permutation_seed(row_seed, permute_dist(x))
  })
}

read_partition_inputs <- function(job) {
  stopifnot(nrow(job) == 1L)
  lapply(job$.files[[1L]], qs2::qs_read, nthreads = 1L)
}

run_partition_job <- function(job, fun, seed, stream = TRUE, inputs = NULL) {
  stopifnot(nrow(job) == 1L)
  if (is.null(inputs)) {
    inputs <- read_partition_inputs(job)
  }
  windows <- sort(unique(unlist(
    lapply(inputs, function(x) {
      if (is.data.frame(x) && "time_id" %in% names(x)) x$time_id else NULL
    }),
    use.names = FALSE
  )))
  evaluate <- function(values) {
    result <- with_permutation_seed(seed, fun(values, seed))
    stopifnot(is.data.frame(result))
    missing <- setdiff(names(job$.key[[1L]]), names(result))
    if (length(missing)) {
      result <- dplyr::bind_cols(
        job$.key[[1L]][rep(1L, nrow(result)), missing, drop = FALSE],
        result
      )
    }
    result
  }
  if (!stream || !length(windows)) {
    return(evaluate(inputs))
  }
  purrr::map(windows, function(time) {
    sliced <- lapply(inputs, function(x) {
      if (is.data.frame(x) && "time_id" %in% names(x)) {
        x[x$time_id == time, , drop = FALSE]
      } else {
        x
      }
    })
    if (
      any(vapply(sliced, function(x) is.data.frame(x) && !nrow(x), logical(1)))
    ) {
      return(NULL)
    }
    evaluate(sliced)
  }) |>
    purrr::list_rbind()
}

run_partition_batch <- function(job, plan, batch, fun) {
  selected <- plan[plan$tar_batch == batch & plan$tar_rep > 0L, ]
  stopifnot(nrow(selected) > 0L)
  inputs <- read_partition_inputs(job)
  purrr::map(seq_len(nrow(selected)), function(i) {
    run_partition_job(job, fun, selected$seed[i], inputs = inputs) |>
      dplyr::mutate(
        tar_batch = selected$tar_batch[i],
        tar_rep = selected$tar_rep[i]
      )
  }) |>
    purrr::list_rbind()
}
