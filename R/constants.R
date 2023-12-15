hypers_alternative <- list(alternative = c("greater", "less"))
hypers_prep_shared <- tidyr::expand_grid(
  resp_trans = c(
    "precise", # use the original response, i.e., 1-4
    "coarse" # collapse 1-2 and 3-4 respectively
    # "compromise" # collapse 3-4 for old and collapse 1-2 for new
  ),
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
hypers_dist_measure <- tibble::tribble(
  ~method, ~call,
  "sm", quote(1 - nomclust::sm(.x)),
  "caylay", quote(calc_dist(.x, Rankcluster::distCayley)),
  "gower", quote(
    .x |>
      mutate(across(everything(), ~ factor(.x, ordered = TRUE))) |>
      proxy::simil(method = "Gower")
  )
)
