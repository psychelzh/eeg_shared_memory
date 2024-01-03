library(targets)
tar_option_set(
  packages = c("tidyverse"),
  format = "qs",
  controller = crew::crew_controller_local(
    name = "local",
    workers = 20
  ),
  memory = "transient"
)
if (Sys.info()["sysname"] == "Windows") {
  future::plan(future.callr::callr)
} else {
  future::plan(future::multicore)
}
tar_source()

# config: check inter-subject similarity ----
if (FALSE) { # we do not need to compare windowed results
  # compare abstract and concrete and subsequent memory effect
  inter_check_window <- tarchetypes::tar_map(
    hypers_rs_window |>
      dplyr::filter(type == "inter"),
    names = c(type, acq, region),
    tar_target(
      rsa_inter_common_trials,
      lapply(
        tar_name_files,
        filter_shared,
        response_shared
      ),
      pattern = map(tar_name_files)
    ),
    tar_target(
      summary_word_cat,
      lapply(
        rsa_inter_common_trials,
        summarise,
        mean_se(fisher_z),
        .by = c(region_id, word_category, window_id)
      ) |>
        list_rbind(),
      pattern = map(rsa_inter_common_trials)
    ),
    tar_target(
      summary_word_mem,
      lapply(
        rsa_inter_common_trials,
        summarise,
        mean_se(fisher_z),
        .by = c(region_id, response_type_shared, window_id)
      ) |>
        list_rbind(),
      pattern = map(rsa_inter_common_trials)
    )
  )
}

# config: predict memory performance ----
targets_pred_perf <- tarchetypes::tar_map(
  hypers_pred_perf,
  tar_target(
    cur_mem_perf,
    filter(mem_perf, .data[["index_name"]] == index_name)
  ),
  tar_target(
    stats_pred_perf,
    extract_stats_pred_perf(avg_rs_group_window, cur_mem_perf)
  ),
  tarchetypes::tar_rep(
    stats_pred_perf_perm,
    extract_stats_pred_perf(
      avg_rs_group_window,
      permutate_behav(cur_mem_perf, "subj_id")
    ),
    batches = num_batches,
    reps = num_reps
  ),
  tar_target(
    clusters_p_pred_perf,
    extract_cluster_p(stats_pred_perf, stats_pred_perf_perm) |>
      add_column(index_name = index_name, .before = 1L)
  )
)

# config: predict shared memory content ----
targets_pred_content <- tarchetypes::tar_map(
  hypers_prep_shared,
  names = c(resp_trans, include),
  tar_target(
    resp_mat,
    events_retrieval |>
      transform_resp() |>
      prepare_resp_mat(include)
  ),
  tarchetypes::tar_map(
    hypers_dist_shared,
    tar_target(
      simil_content,
      calc_dist_resp_mat(resp_mat, method = method),
      deployment = "main"
    ),
    tarchetypes::tar_map(
      tibble::tribble(
        ~mantel, ~covariate,
        "mantel", NULL,
        "partial", quote(dist_mem_perf)
      ),
      names = mantel,
      tar_target(
        stats_pred_content,
        extract_stats_pred_content(
          avg_rs_inter_trial,
          simil_content,
          covariate = covariate,
          keep_perms = TRUE
        )
      )
    )
  )
)

# main targets definition ----
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
    tarchetypes::tar_files_input(
      tar_name_path,
      config_files_rs(type, acq),
      batches = batches_file
    ),
    hypers_rs
  ),
  tarchetypes::tar_map(
    dplyr::filter(hypers_rs, acq != "whole"),
    names = c(type, acq),
    tar_target(
      avg_rs,
      lapply(
        tar_name_path,
        average_rs_trials,
        scalar_rs = type == "group"
      ) |>
        list_rbind(),
      pattern = map(tar_name_path)
    )
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
  tar_target(
    rsa_inter_avg_by_category,
    rsa_inter_common_trials |>
      summarise(
        mean_fisher_z = mean(fisher_z, na.rm = TRUE),
        .by = c(region_id, subj_id_col, subj_id_row, word_category)
      )
  ),
  tar_target(
    rsa_inter_avg_by_resp,
    rsa_inter_common_trials |>
      summarise(
        mean_fisher_z = mean(fisher_z, na.rm = TRUE),
        .by = c(region_id, subj_id_col, subj_id_row, response_type_shared)
      )
  ),
  tar_target(
    stats_rsa_inter_by_resp,
    extract_stats_sme(rsa_inter_avg_by_resp)
  ),
  if (FALSE) list(
    inter_check_window,
    tarchetypes::tar_combine(
      summary_word_cat_rsa_inter_common_trials_window,
      inter_check_window$summary_word_cat
    ),
    tarchetypes::tar_combine(
      summary_word_mem_rsa_inter_common_trials_window,
      inter_check_window$summary_word_mem
    )
  ),
  # predict memory performance ----
  tar_target(mem_perf, calc_mem_perf(events_retrieval)),
  tar_target(dist_mem_perf, calc_dist_mem_perf(mem_perf)),
  targets_pred_perf,
  tarchetypes::tar_combine(
    stats_pred_perf,
    targets_pred_perf$stats_pred_perf
  ),
  tarchetypes::tar_combine(
    clusters_p_pred_perf,
    targets_pred_perf$clusters_p_pred_perf
  ),
  # predict shared memory content ----
  targets_pred_content,
  tar_combine_with_meta(
    stats_pred_content,
    select_list(
      targets_pred_content,
      starts_with("stats_pred_content")
    ),
    cols_targets = c("mantel", "method", "resp_trans", "include"),
    prefix = "stats_pred_content"
  ),
  tar_target(
    file_stats_pred_content_real,
    "data/pred_content/res_par-mantel_isc_rps_trial_avg_smc_gower.csv",
    format = "file"
  ),
  tarchetypes::tar_group_by(
    stats_pred_content_real,
    read_csv(file_stats_pred_content_real, show_col_types = FALSE) |>
      rename(statistic.r = statistic_r),
    mantel_type, method, include
  ),
  tar_target(
    file_stats_pred_content_perm,
    "data/pred_content/res_par-mantel_isc_rps_trial-avg_smc_rand1000_gower.csv",
    format = "file"
  ),
  tarchetypes::tar_group_by(
    stats_pred_content_perm,
    read_csv(file_stats_pred_content_perm, show_col_types = FALSE),
    mantel_type, method, include
  ),
  tar_target(
    clusters_p_pred_content,
    extract_cluster_p(
      stats_pred_content_real,
      stats_pred_content_perm,
      cols_group = c(method, include, mantel_type),
      cols_perm = perm_id,
      col_statistic = statistic.r
    ),
    pattern = map(stats_pred_content_real, stats_pred_content_perm)
  ),
  # render website ----
  tarchetypes::tar_quarto(website)
)
