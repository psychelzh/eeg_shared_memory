tar_mantel <- function(
  name,
  data_whole,
  data_dynamic,
  ydis,
  zdis = NULL,
  reps = 10L,
  batches = 100L,
  permutations = 9999L
) {
  tar_whole <- function(name, data_expr) {
    name_data <- paste0("data_", name, "_whole")
    name_stats <- paste0("stats_", name, "_whole")
    list(
      eval(bquote(tar_partition_target(
        .(name_data),
        .(substitute(data_expr))
      ))),
      tar_target_raw(
        name_stats,
        bquote(extract_stats_mantel(.(as.symbol(name_data))))
      )
    )
  }
  if (is.null(substitute(zdis))) {
    data_whole_expr <- substitute(calc_mantel(
      data_whole,
      ydis,
      permutations = permutations
    ))
    data_dynamic_expr <- substitute(calc_mantel(
      data_dynamic,
      ydis,
      permutations = permutations
    ))
    data_dynamic_perm_expr <- substitute(
      calc_mantel(data_dynamic, permute_dist(ydis), permutations = permutations)
    )
  } else {
    data_whole_expr <- substitute(calc_mantel_partial(
      data_whole,
      ydis,
      zdis,
      permutations = permutations
    ))
    data_dynamic_expr <- substitute(calc_mantel_partial(
      data_dynamic,
      ydis,
      zdis,
      permutations = permutations
    ))
    data_dynamic_perm_expr <- substitute(
      calc_mantel_partial(
        data_dynamic,
        permute_dist(ydis),
        zdis,
        permutations = permutations
      )
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
          reps = .(reps),
          batches = .(batches)
        )
      )
    )
  )
}
