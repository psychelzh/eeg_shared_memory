# paths and files ----
config_path_dataset <- function(type, acq, region) {
  fs::path(
    "data",
    sprintf("type-%s_acq-%s", type, acq),
    sprintf("region-%s", region)
  )
}

config_path_file <- function(type, acq) {
  fs::path(
    "data",
    sprintf(
      "type-%s_acq-%s_rs.parquet",
      type, acq
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
