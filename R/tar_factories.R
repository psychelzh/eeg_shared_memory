tar_cluster_permutation <- function(
  name,
  stats_expr,
  stats_perm_expr,
  data_expr = NULL,
  data_perm_expr = NULL,
  clusters_stats_expr = NULL,
  stats_name = NULL,
  data_name = NULL,
  clusters_stats_name = NULL
) {
  if (!missing(name)) {
    stats_name <- paste0("stats_", name)
    data_name <- paste0("data_", name)
    clusters_stats_name <- paste0("clusters_stats_", name)
  }
  stats_name_permuted <- paste0(stats_name, "_permuted")
  data_name_permuted <- paste0(data_name, "_permuted")
  if (!is.null(substitute(data_expr))) {
    stats_expr <- targets::tar_tidy_eval(
      as.expression(substitute(stats_expr)),
      envir = list(.x = rlang::sym(data_name)),
      tidy_eval = TRUE
    )
    stats_perm_expr <- targets::tar_tidy_eval(
      as.expression(substitute(stats_perm_expr)),
      envir = list(.x = rlang::sym(data_name_permuted)),
      tidy_eval = TRUE
    )
  }
  if (is.null(substitute(clusters_stats_expr))) {
    clusters_stats_expr <- rlang::call2(
      "calc_clusters_stats",
      rlang::sym(stats_name),
      rlang::sym(stats_name_permuted)
    )
  } else {
    clusters_stats_expr <- targets::tar_tidy_eval(
      as.expression(substitute(clusters_stats_expr)),
      envir = list(
        .x = rlang::sym(stats_name),
        .y = rlang::sym(stats_name_permuted)
      ),
      tidy_eval = TRUE
    )
  }
  list(
    targets::tar_target_raw(stats_name, substitute(stats_expr)),
    if (!is.null(substitute(data_expr))) {
      list(
        targets::tar_target_raw(data_name, substitute(data_expr)),
        tarchetypes::tar_rep_raw(
          data_name_permuted,
          substitute(data_perm_expr),
          reps = 10,
          batches = 100
        ),
        tarchetypes::tar_rep2_raw(
          stats_name_permuted,
          substitute(stats_perm_expr),
          data_name_permuted
        )
      )
    } else {
      tarchetypes::tar_rep_raw(
        stats_name_permuted,
        substitute(stats_perm_expr),
        reps = 10,
        batches = 100
      )
    },
    targets::tar_target_raw(clusters_stats_name, clusters_stats_expr)
  )
}

tar_mantel <- function(name, data_whole, data_dynamic, ydis, zdis = NULL) {
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
          stats_perm_expr = extract_stats_mantel(!!.x)
        )
      )
    )
  )
}
