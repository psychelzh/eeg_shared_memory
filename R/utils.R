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

select_list <- function(.l, ...) {
  pos <- tidyselect::eval_select(
    rlang::expr(c(...)),
    .l
  )
  rlang::set_names(.l[pos], names(pos))
}
