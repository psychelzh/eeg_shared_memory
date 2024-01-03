# paths and files ----
config_files_rs <- function(type, acq) {
  switch(acq,
         window = fs::path(
           "data",
           sprintf("type-%s_acq-%s", type, acq)
         ) |>
           fs::dir_ls(
             recurse = TRUE,
             type = "file"
           ),
         trial = ,
         whole = fs::path(
           "data",
           sprintf(
             "type-%s_acq-%s_rs.parquet",
             type, acq
           )
         ),
         stop("Unknown parameter")
  )
}

config_files_pred_perf_rps <- function(index_name_sjt) {
  fs::path(
    "data/representational_space",
    sprintf(
      "corr_indiv2grp_rps_spc_%s_rightside.mat",
      index_name_sjt
    )
  )
}

config_files_pred_content_rps <- function(method, type = c("real", "perm")) {
  fs::path(
    "data/representational_space",
    sprintf(
      "res_par-mantel_isc_rps_spc_memory_content%s_%s.csv",
      if_else(type == "real", "", "_rand1000"), method
    )
  )
}

# programming ----
select_list <- function(.l, ...) {
  pos <- tidyselect::eval_select(
    rlang::expr(c(...)),
    .l
  )
  rlang::set_names(.l[pos], names(pos))
}

# distance ----
as_dist_vec <- function(vec, ..., size = NULL, diag = FALSE, upper = FALSE) {
  if (is.null(size)) {
    size <- 0.5 + sqrt(0.25 + 2 * length(vec))
  }
  stopifnot(all.equal(size, as.integer(size)))
  structure(
    vec,
    class = "dist",
    Size = size,
    Diag = diag,
    Upper = upper,
    ...
  )
}

stretch <- function(mat, name_value = "val") {
  d <- as.dist(mat)
  combn(names(d), 2, simplify = FALSE) |>
    do.call(rbind, args = _) |>
    as_tibble(.name_repair = ~ c("x", "y")) |>
    mutate("{name_value}" := as.vector(d))
}
