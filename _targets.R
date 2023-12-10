library(targets)
tar_option_set(
  packages = c("tidyverse"),
  format = "qs",
  controller = crew::crew_controller_local(
    name = "local",
    workers = 8
  )
)
future::plan(future.callr::callr)
tar_source()

config_window_rs <- tidyr::expand_grid(
  type = c("inter", "group"),
  acq = c("window"),
  region = paste0("region", 1:6)
)

group_pred_perf <- tarchetypes::tar_map(
  config_window_rs |>
    dplyr::filter(type == "group"),
  tar_target(
    path_dataset,
    config_path_dataset(type, acq, region),
    format = "file"
  ),
  tar_target(
    stats,
    extract_stats_group(path_dataset, mem_perf)
  ),
  tarchetypes::tar_rep(
    stats_perm,
    extract_stats_group(
      path_dataset,
      permutate_behav(mem_perf, "subj_id")
    ),
    batches = 100,
    reps = 10
  )
)

list(
  tarchetypes::tar_file_read(
    events_encoding,
    "data/group_task-wordencoding_events.csv",
    read = readr::read_csv(!!.x, show_col_types = FALSE)
  ),
  tarchetypes::tar_file_read(
    events_retrieval,
    "data/group_task-wordretrieval_events.csv",
    read = readr::read_csv(!!.x, show_col_types = FALSE)
  ),
  tar_target(subj_pair_filter, prepare_subj_pair_common(events_retrieval)),
  tar_target(mem_perf, calc_mem_perf(events_retrieval)),
  group_pred_perf
)
