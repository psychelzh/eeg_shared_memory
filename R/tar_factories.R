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
  pattern = NULL,
  reps = 1000
) {
  if (!missing(name)) {
    stats_name <- paste0("stats_", name)
    data_name <- paste0("data_", name)
    clusters_stats_name <- paste0("clusters_stats_", name)
  }
  stats_name_permuted <- paste0(stats_name, "_permuted")
  data_name_permuted <- paste0(data_name, "_permuted")
  if (!is.null(substitute(data_expr))) {
    stats_expr <- tar_tidy_eval(
      as.expression(substitute(stats_expr)),
      envir = list(.x = rlang::sym(data_name)),
      tidy_eval = TRUE
    )
    stats_perm_expr <- tar_tidy_eval(
      as.expression(substitute(stats_perm_expr)),
      envir = list(.x = rlang::sym(data_name_permuted)),
      tidy_eval = TRUE
    )
    if (!is.null(substitute(pattern))) {
      pattern_data <- substitute(pattern)
      pattern_stats <- rlang::call2("map", rlang::sym(data_name))
      pattern_stats_permuted <- rlang::call2(
        "map",
        rlang::sym(data_name_permuted)
      )
    } else {
      pattern_data <- pattern_stats <- pattern_stats_permuted <- NULL
    }
  } else {
    pattern_stats <- pattern_stats_permuted <- substitute(pattern)
  }
  if (is.null(substitute(clusters_stats_expr))) {
    clusters_stats_expr <- rlang::call2(
      "calc_clusters_stats",
      rlang::sym(stats_name),
      rlang::sym(stats_name_permuted)
    )
  } else {
    clusters_stats_expr <- tar_tidy_eval(
      as.expression(substitute(clusters_stats_expr)),
      envir = list(
        .x = rlang::sym(stats_name),
        .y = rlang::sym(stats_name_permuted)
      ),
      tidy_eval = TRUE
    )
  }
  list(
    tar_target_raw(stats_name, substitute(stats_expr), pattern = pattern_stats),
    if (!is.null(substitute(data_expr))) {
      list(
        tar_target_raw(
          data_name,
          substitute(data_expr),
          pattern = pattern_data,
          iteration = "list"
        ),
        tar_target_raw(
          data_name_permuted,
          bquote(run_rep_df(.(substitute(data_perm_expr)), .(reps))),
          pattern = pattern_data,
          iteration = "list"
        ),
        tar_target_raw(
          stats_name_permuted,
          substitute(stats_perm_expr),
          pattern = pattern_stats_permuted
        )
      )
    } else {
      tarchetypes::tar_rep_raw(
        stats_name_permuted,
        bquote(run_rep_df(.(substitute(stats_perm_expr)), .(reps))),
        pattern = pattern_stats_permuted
      )
    },
    tar_target_raw(clusters_stats_name, clusters_stats_expr)
  )
}

tar_mantel <- function(
  name,
  data_whole,
  data_dynamic,
  ydis,
  zdis = NULL,
  pattern = NULL,
  reps = 1000
) {
  tar_whole <- function(name, data_expr) {
    name_data <- paste0("data_", name, "_whole")
    name_stats <- paste0("stats_", name, "_whole")
    list(
      tar_target_raw(name_data, substitute(data_expr)),
      tar_target_raw(
        name_stats,
        bquote(extract_stats_mantel(.(as.symbol(name_data))))
      )
    )
  }
  if (is.null(substitute(zdis))) {
    data_whole_expr <- substitute(calc_mantel(data_whole, ydis))
    data_dynamic_expr <- substitute(calc_mantel(data_dynamic, ydis))
    data_dynamic_perm_expr <- substitute(
      calc_mantel(data_dynamic, permute_dist(ydis))
    )
  } else {
    data_whole_expr <- substitute(calc_mantel_partial(data_whole, ydis, zdis))
    data_dynamic_expr <- substitute(calc_mantel_partial(
      data_dynamic,
      ydis,
      zdis
    ))
    data_dynamic_perm_expr <- substitute(
      calc_mantel_partial(data_dynamic, permute_dist(ydis), zdis)
    )
  }
  eval(
    bquote(
      list(
        tar_whole(name, .(data_whole_expr)),
        tar_cluster_permutation(
          paste0(name, "_dynamic"),
          data_expr = .(data_dynamic_expr),
          data_perm_expr = .(data_dynamic_perm_expr),
          stats_expr = extract_stats_mantel(!!.x),
          stats_perm_expr = extract_stats_mantel(!!.x),
          pattern = .(substitute(pattern)),
          reps = .(reps)
        )
      )
    )
  )
}

run_rep_df <- function(expr, n) {
  replicate(n, expr, simplify = FALSE) |>
    list_rbind(names_to = ".rep")
}
