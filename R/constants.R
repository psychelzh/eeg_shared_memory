# permutations number: 1000 divided into 100 batches of 10 reps
num_batches <- 100
num_reps <- 10

hypers_pred_perf <- tibble::tribble(
  ~index_name, ~index_name_sjt,
  "remember", "rem_dp",
  "knowadj", "know_dp_adj",
  "avg_rk", "mean_rem_know",
  "avg_score", "sum_mean"
)

hypers_prep_shared <- tidyr::expand_grid(
  resp_trans = "precise",
  # c(
  #   "precise", # use the original response, i.e., 1-4
  #   "coarse", # collapse 1-2 and 3-4 respectively
  #   "compromise" # collapse 3-4 for old and collapse 1-2 for new
  # ),
  include = c("old", "all")
) |>
  dplyr::mutate(
    transform_resp = rlang::syms(
      sprintf("transform_resp_%s", resp_trans)
    ),
    tar_resp_mat = rlang::syms(
      sprintf("resp_mat_%s_%s", resp_trans, include)
    )
  )
hypers_dist_shared <- tibble::tibble(method = c("sm", "gower"))

hypers_rs <- tidyr::expand_grid(
  type = c("inter", "group"),
  acq = c("trial", "whole", "window")
) |>
  dplyr::mutate(
    tar_name_path = rlang::syms(
      sprintf("file_rs_%s_%s", type, acq)
    ),
    batches_file = dplyr::if_else(acq == "window", 50, 1)
  )

labels_acq <- c(
  "trial" = "Averaged Trial-level",
  "whole" = "Concatenated Time-series"
)
labels_index_name <- c(
  "knowadj" = "Familiarity",
  "remember" = "Recollection",
  "avg_rk" = "Average F/R",
  "avg_score" = "Memory Score"
)
labels_include <- c(
  "all" = "Both Items",
  "old" = "Old Items Only"
)
labels_method <- c(
  "gower" = "Gower (ordinal)",
  "sm" = "Simple Match (nomial)"
)
