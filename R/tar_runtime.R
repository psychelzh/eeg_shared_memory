# Resolve a schedulable memory resource before any SGE jobs are submitted.
select_sge_memory_config <- function(
  lines,
  resource = "",
  min_host_gigabytes = NULL
) {
  lines <- lines[!grepl("^\\s*(#|$)", lines)]
  fields <- strsplit(trimws(lines), "[[:space:]]+")
  fields <- fields[lengths(fields) >= 6L]
  eligible <- vapply(
    fields,
    function(x) {
      x[3L] == "MEMORY" && x[4L] == "<=" && x[5L] %in% c("YES", "FORCED")
    },
    logical(1)
  )
  available <- vapply(fields[eligible], `[[`, character(1), 1L)
  candidates <- if (nzchar(resource)) resource else c("m_mem_free", "mem_free")
  selected <- candidates[candidates %in% available]
  if (!length(selected)) {
    stop(
      "No matching SGE memory resource. Inspect qconf -sc and set ",
      "EEG_SGE_MEMORY_RESOURCE to a requestable MEMORY resource with relop <=. ",
      "Available: ",
      paste(available, collapse = ", "),
      ". No SGE jobs have been submitted."
    )
  }
  selected <- selected[1L]
  consumable <- fields[eligible][[match(selected, available)]][6L] %in%
    c("YES", "JOB")
  if (is.null(min_host_gigabytes)) {
    min_host_gigabytes <- if (consumable) 0 else 400
  }
  if (
    length(min_host_gigabytes) != 1L ||
      !is.finite(min_host_gigabytes) ||
      min_host_gigabytes < 0
  ) {
    stop("EEG_SGE_MIN_HOST_GB must be a nonnegative number.")
  }
  if (min_host_gigabytes > 0 && !"mem_total" %in% available) {
    stop(
      "The minimum-host-memory filter requires a requestable mem_total resource."
    )
  }
  list(
    resource = selected,
    consumable = consumable,
    min_host_gigabytes = min_host_gigabytes
  )
}

read_sge_memory_config <- function() {
  qconf <- Sys.which("qconf")
  if (!nzchar(qconf)) {
    stop("SGE configuration requires qconf on PATH.")
  }
  lines <- system2(qconf, "-sc", stdout = TRUE, stderr = TRUE)
  status <- attr(lines, "status")
  if (!is.null(status) && status != 0L) {
    stop("qconf -sc failed: ", paste(lines, collapse = "\n"))
  }
  min_host <- Sys.getenv("EEG_SGE_MIN_HOST_GB", "")
  config <- select_sge_memory_config(
    lines,
    Sys.getenv("EEG_SGE_MEMORY_RESOURCE", ""),
    if (nzchar(min_host)) suppressWarnings(as.numeric(min_host)) else NULL
  )
  if (!config$consumable) {
    message(
      "SGE ",
      config$resource,
      " is not consumable: requests are free-memory ",
      "thresholds, not reservations. Minimum host memory: ",
      config$min_host_gigabytes,
      " GB."
    )
  }
  config
}

sge_memory_directive <- function(resource, gigabytes) {
  stopifnot(
    length(resource) == 1L,
    grepl("^[A-Za-z][A-Za-z0-9_]*$", resource),
    length(gigabytes) == 1L,
    is.finite(gigabytes),
    gigabytes > 0
  )
  paste0("#$ -l ", resource, "=", gigabytes, "G")
}

# Configure execution without adding scheduler variables to the target globals.
configure_targets_runtime <- function(
  on_sge = tolower(Sys.info()["nodename"]) %in% c("shadow", "hippocampus"),
  memory_config = NULL
) {
  # Local default: i7-12700F (12 cores / 20 threads), 64 GiB RAM.
  # Prefer stability after a full-run worker startup failure at 14 workers.
  # Keep headroom for preparation, scheduling and interactive applications.
  worker_count_variable <- if (on_sge) "EEG_SGE_WORKERS" else "EEG_WORKERS"
  configured_workers <- Sys.getenv(
    worker_count_variable,
    if (on_sge) "32" else "12"
  )
  compute_workers <- suppressWarnings(as.integer(configured_workers))
  if (
    is.na(compute_workers) ||
      compute_workers < 1L ||
      as.character(compute_workers) != configured_workers
  ) {
    stop(worker_count_variable, " must be a positive integer.")
  }
  # Skip slow system-library sandbox setup; retain project renv activation.
  # Set limits before startup, including explicitly in remote SGE job scripts.
  worker_environment <- c(
    RENV_CONFIG_SANDBOX_ENABLED = "FALSE",
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    RCPP_PARALLEL_NUM_THREADS = "1"
  )
  do.call(Sys.setenv, as.list(worker_environment))
  if (on_sge && is.null(memory_config)) {
    memory_config <- read_sge_memory_config()
  }
  # File preparation remains serial in its own pool.
  controller_pool <- function(name, workers) {
    if (on_sge) {
      log_directory <- file.path(getwd(), "logs", "crew", "sge", name)
      script_directory <- file.path(log_directory, "scripts")
      dir.create(script_directory, recursive = TRUE, showWarnings = FALSE)
      crew.cluster::crew_controller_sge(
        name = name,
        workers = workers,
        seconds_idle = 300,
        options_cluster = crew.cluster::crew_options_sge(
          script_directory = script_directory,
          # A serial SGE job uses one slot; do not assume a parallel environment.
          script_lines = c(
            sge_memory_directive(
              memory_config$resource,
              if (name == "heavy") 24 else 6
            ),
            if (memory_config$min_host_gigabytes > 0) {
              sge_memory_directive(
                "mem_total",
                memory_config$min_host_gigabytes
              )
            },
            paste0(
              "export ",
              names(worker_environment),
              "=",
              shQuote(worker_environment, type = "sh")
            )
          ),
          log_output = paste0(
            normalizePath(log_directory, winslash = "/", mustWork = TRUE),
            "/"
          )
        )
      )
    } else {
      crew::crew_controller_local(
        name = name,
        workers = workers,
        seconds_launch = 120,
        options_local = crew::crew_options_local(
          log_directory = file.path(getwd(), "logs", "crew", name)
        )
      )
    }
  }
  controller <- crew::crew_controller_group(
    controller_pool("light", compute_workers),
    controller_pool("heavy", 1)
  )
  targets::tar_option_set(
    packages = c("tidyverse"),
    format = "qs",
    resources = targets::tar_resources(
      crew = targets::tar_resources_crew(controller = "light")
    ),
    controller = controller
  )
  invisible(controller)
}
