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
