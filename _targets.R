library(targets)
tar_option_set(
  packages = c("tidyverse"),
  format = "qs",
  controller = crew::crew_controller_local(
    name = "local",
    workers = 20
  ),
  memory = "transient",
  garbage_collection = TRUE
)
future::plan(future.callr::callr)
tar_source()

config_nonwin_rs <- tidyr::expand_grid(
  type = c("inter", "group"),
  acq = c("trial", "whole")
) |>
  dplyr::mutate(
    tar_name_file = rlang::syms(
      sprintf("file_rs_%s_%s", type, acq)
    )
  )
config_window_rs <- tidyr::expand_grid(
  type = c("inter", "group"),
  acq = c("window"),
  region = paste0("region", 1:6)
) |>
  dplyr::mutate(
    tar_name_path = rlang::syms(
      sprintf("path_dataset_%s_%s_%s", type, acq, region)
    )
  )

group_pred_perf <- tarchetypes::tar_map(
  config_window_rs |>
    dplyr::filter(type == "group"),
  names = c(type, acq, region),
  tar_target(
    stats,
    extract_stats_group(tar_name_path, mem_perf)
  ),
  tarchetypes::tar_rep(
    stats_perm,
    extract_stats_group(
      tar_name_path,
      permutate_behav(mem_perf, "subj_id")
    ),
    batches = 100,
    reps = 10
  )
)

inter_check_window <- tarchetypes::tar_map(
  config_window_rs |>
    dplyr::filter(type == "inter"),
  names = c(type, acq, region),
  tar_target(
    path_chunks,
    fs::dir_ls(tar_name_path, type = "file", recurse = TRUE) |>
      as_tibble_col()
  ),
  tar_target(
    rsa_inter_common_trials,
    lapply(
      unlist(path_chunks),
      filter_inter_rs_by_trial,
      events_encoding,
      subj_pair_filter
    ),
    pattern = map(path_chunks)
  ),
  tar_target(
    summary_word_cat_branches,
    rsa_inter_common_trials |>
      lapply(
        summarise,
        n = n(),
        sum_fisher_z = sum(fisher_z),
        .by = c(region_id, word_category, window_id)
      ) |>
      list_rbind(),
    pattern = map(rsa_inter_common_trials)
  ),
  tar_target(
    summary_word_cat,
    summary_word_cat_branches |>
      summarise(
        mean_fisher_z = sum(sum_fisher_z) / sum(n),
        .by = c(region_id, word_category, window_id)
      )
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
  tarchetypes::tar_eval(
    tar_target(
      tar_name_file,
      config_path_file(type, acq),
      format = "file"
    ),
    values = config_nonwin_rs
  ),
  tarchetypes::tar_eval(
    tar_target(
      tar_name_path,
      config_path_dataset(type, acq, region),
      format = "file_fast"
    ),
    values = config_window_rs
  ),
  tar_target(
    subj_pair_filter,
    prepare_subj_pair_common(events_retrieval)
  ),
  tar_target(
    rsa_inter_common_trials,
    filter_inter_rs_by_trial(
      file_rs_inter_trial,
      events_encoding,
      subj_pair_filter
    )
  ),
  tar_target(mem_perf, calc_mem_perf(events_retrieval)),
  group_pred_perf,
  inter_check_window,
  tarchetypes::tar_combine(
    summary_word_cat_rsa_inter_common_trials_window,
    inter_check_window$summary_word_cat
  ),
  tarchetypes::tar_quarto(website)
)
