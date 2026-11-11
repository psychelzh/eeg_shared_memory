# Native partitions are persisted as individual targets branches. Catalogs and
# jobs contain keys only; the branch payload arrives through targets retrieval.
native_partition_keys <- function(data, by = "cca_id") {
  stopifnot(is.data.frame(data), length(by) > 0L, all(by %in% names(data)))
  keys <- data |>
    dplyr::ungroup() |>
    tibble::as_tibble() |>
    dplyr::distinct(dplyr::across(dplyr::all_of(by))) |>
    dplyr::arrange(dplyr::across(dplyr::everything()))
  stopifnot(nrow(keys) > 0L, !anyNA(keys), !anyDuplicated(keys))
  keys
}

native_partition_piece <- function(data, key) {
  stopifnot(nrow(key) == 1L, !anyNA(key), all(names(key) %in% names(data)))
  piece <- dplyr::semi_join(data, key, by = names(key))
  stopifnot(nrow(piece) > 0L)
  list(key = key, data = piece)
}

native_eeg_partition_piece <- function(file, channel_regions, key) {
  stopifnot(nrow(key) == 1L, key$region_id %in% seq_len(6L))
  data <- qs2::qs_read(file, nthreads = 1L)
  channels <- channel_regions[[paste0("region", key$region_id)]] != 0
  list(key = key, data = data[channels, , , , drop = FALSE])
}

make_native_partition_jobs <- function(catalogs) {
  if (!length(catalogs)) {
    return(tibble::tibble(.key = list(tibble::tibble())))
  }
  reference <- catalogs[[1L]]
  stopifnot(nrow(reference) > 0L, !anyNA(reference), !anyDuplicated(reference))
  # map() zips branch positions. Reject both missing keys and a changed order.
  for (catalog in catalogs) {
    if (!identical(reference, catalog)) {
      stop("Native partition sources must have identical ordered key catalogs.")
    }
  }
  tibble::tibble(
    .key = lapply(seq_len(nrow(reference)), function(i) {
      reference[i, , drop = FALSE]
    })
  )
}

native_partition_inputs <- function(job, pieces) {
  stopifnot(nrow(job) == 1L)
  lapply(pieces, function(piece) {
    if (!identical(job$.key[[1L]], piece$key)) {
      stop("Native partition key does not match its computation job.")
    }
    piece$data
  })
}

# Legacy file helpers remain for old-cache reading and reference comparisons.
# The production target factories below use native objects, never these bundles.
partition_relative_path <- function(path) {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  absolute <- normalizePath(path, winslash = "/", mustWork = TRUE)
  vapply(
    absolute,
    function(file) {
      prefix <- paste0(root, "/")
      inside <- if (.Platform$OS.type == "windows") {
        startsWith(tolower(file), tolower(prefix))
      } else {
        startsWith(file, prefix)
      }
      # External stores remain absolute. In particular, a C: temporary store
      # cannot be represented as a relative path from a D: project on Windows.
      if (inside) as.character(fs::path_rel(file, start = root)) else file
    },
    character(1),
    USE.NAMES = FALSE
  )
}

# Older bundles and jobs recorded a machine's absolute project root. Prefer the
# corresponding local store when reading such a copied cache. New outputs below
# never record that root. Paths outside the default store are left unchanged.
partition_file_path <- function(path) {
  path <- gsub("\\\\", "/", path)
  vapply(
    path,
    function(file) {
      absolute <- grepl("^([A-Za-z]:/|/)", file)
      marker <- regexpr("/_targets/", file, fixed = TRUE)[1L]
      if (absolute && marker > 0L) {
        local <- substring(file, marker + 1L)
        if (file.exists(local)) return(local)
      }
      file
    },
    character(1),
    USE.NAMES = FALSE
  )
}

write_partition_bundle <- function(keys, get_piece, directory) {
  stopifnot(nrow(keys) > 0L, !anyNA(keys), !anyDuplicated(keys))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  staging <- tempfile("staging-", tmpdir = directory)
  dir.create(staging)
  # Cleanup is restricted to this function's newly created staging directory.
  on.exit(
    {
      stopifnot(identical(
        dirname(normalizePath(staging, winslash = "/", mustWork = TRUE)),
        directory
      ))
      unlink(staging, recursive = TRUE)
    },
    add = TRUE
  )
  files <- vapply(
    seq_len(nrow(keys)),
    function(i) {
      file <- file.path(staging, paste0("part-", i, ".qs2"))
      qs2::qs_save(get_piece(i), file, nthreads = 1L)
      file
    },
    character(1)
  )
  hashes <- vapply(
    files,
    digest::digest,
    character(1),
    file = TRUE,
    algo = "xxhash64"
  )
  signature <- digest::digest(
    list(version = 1L, keys = keys, hashes = unname(hashes)),
    algo = "xxhash64",
    serializeVersion = 2L
  )
  generation <- file.path(directory, paste0("generation-", signature))
  dir.create(generation, showWarnings = FALSE)
  manifest <- file.path(staging, "manifest.rds")
  # Manifest contents are independent of the project root and generation path.
  saveRDS(dplyr::mutate(keys, .file = basename(files)), manifest, version = 2L)
  staged <- c(manifest, files)
  outputs <- file.path(generation, basename(staged))
  for (i in seq_along(staged)) {
    unchanged <- file.exists(outputs[i]) &&
      identical(
        digest::digest(file = outputs[i], algo = "xxhash64"),
        digest::digest(file = staged[i], algo = "xxhash64")
      )
    if (!unchanged && !file.copy(staged[i], outputs[i], overwrite = TRUE)) {
      stop("Could not write partition file: ", outputs[i])
    }
  }
  partition_relative_path(outputs)
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
  manifests <- lapply(bundles, function(paths) {
    manifest <- partition_file_path(paths[1L])
    part <- readRDS(manifest)
    # Both legacy absolute entries and new filename-only entries refer to pieces
    # beside this manifest, not to the root of the machine that produced it.
    files <- basename(gsub("\\\\", "/", part$.file))
    part$.file <- partition_relative_path(file.path(dirname(manifest), files))
    part
  })
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
  lapply(job$.files[[1L]], function(file) {
    qs2::qs_read(partition_file_path(file), nthreads = 1L)
  })
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

run_partition_batch <- function(job, plan, batch, fun, inputs = NULL) {
  selected <- plan[plan$tar_batch == batch & plan$tar_rep > 0L, ]
  stopifnot(nrow(selected) > 0L)
  if (is.null(inputs)) {
    inputs <- read_partition_inputs(job)
  }
  purrr::map(seq_len(nrow(selected)), function(i) {
    run_partition_job(job, fun, selected$seed[i], inputs = inputs) |>
      dplyr::mutate(
        tar_batch = selected$tar_batch[i],
        tar_rep = selected$tar_rep[i]
      )
  }) |>
    purrr::list_rbind()
}
