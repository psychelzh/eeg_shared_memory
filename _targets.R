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
    ),
    tar_path_chunks = rlang::syms(
      sprintf("path_chunks_%s_%s_%s", type, acq, region)
    )
  )

# compare abstract and concrete and subsequent memory effect
inter_check_window <- tarchetypes::tar_map(
  config_window_rs |>
    dplyr::filter(type == "inter"),
  names = c(type, acq, region),
  tar_target(
    rsa_inter_common_trials,
    lapply(
      unlist(tar_path_chunks),
      filter_shared,
      response_shared
    ),
    pattern = map(tar_path_chunks)
  ),
  tar_summary_with_branches(
    summary_word_cat,
    rsa_inter_common_trials,
    .by = c(region_id, word_category, window_id)
  ),
  tar_summary_with_branches(
    summary_word_mem,
    rsa_inter_common_trials,
    .by = c(region_id, response_type_shared, window_id)
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
  ),
  tarchetypes::tar_map(
    hypers_alternative,
    tar_target(
      clusters,
      extract_cluster_stats(stats, p.value, statistic, alternative)
    ),
    tar_target(
      clusters_perm,
      extract_cluster_stats(stats_perm, p.value, statistic, alternative)
    ),
    tar_target(
      clusters_p,
      clusters |>
        left_join(
          clusters_perm,
          by = c("region_id", "mem_type"),
          suffix = c("", "_perm")
        ) |>
        summarise(
          p_perm = mean(sum_t_perm > sum_t),
          .by = c(region_id, mem_type, start, end, sum_t)
        )
    )
  )
)

shared_content <- tarchetypes::tar_map(
  hypers_prep_shared,
  names = c(resp_trans, include),
  tar_target(
    resp_mat,
    events_retrieval |>
      transform_resp() |>
      prepare_resp_mat(include)
  ),
  tarchetypes::tar_map(
    hypers_dist_measure,
    names = method,
    tar_target(
      simil,
      eval(substitute(call), envir = list(.x = resp_mat))
    )
  )
)

list(
  # prepare files and paths ----
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
  tarchetypes::tar_eval(
    tar_target(
      tar_path_chunks,
      fs::dir_ls(tar_name_path, type = "file", recurse = TRUE) |>
        split(1:50)
    ),
    values = config_window_rs
  ),
  # check inter-subject similarity ----
  tar_target(
    response_shared,
    extract_response_shared(events_encoding, events_retrieval)
  ),
  tar_target(
    rsa_inter_common_trials,
    filter_shared(file_rs_inter_trial, response_shared)
  ),
  inter_check_window,
  tarchetypes::tar_combine(
    summary_word_cat_rsa_inter_common_trials_window,
    inter_check_window$summary_word_cat
  ),
  tarchetypes::tar_combine(
    summary_word_mem_rsa_inter_common_trials_window,
    inter_check_window$summary_word_mem
  ),
  # predict memory performance ----
  tar_target(mem_perf, calc_mem_perf(events_retrieval)),
  group_pred_perf,
  tarchetypes::tar_combine(
    stats_group_window,
    group_pred_perf$stats
  ),
  tarchetypes::tar_combine(
    clusters_p_greater_group_window,
    group_pred_perf$clusters_p_greater
  ),
  tarchetypes::tar_combine(
    clusters_p_less_group_window,
    group_pred_perf$clusters_p_less
  ),
  # predict shared memory content ----
  shared_content,
  tar_combine_with_meta(
    simil,
    select_list(shared_content, starts_with("simil")),
    cols_targets = c("method", "resp_trans", "include"),
    fun_pre = ~ tibble(mat = list(.x))
  ),
  # render website ----
  tarchetypes::tar_quarto(website)
)
