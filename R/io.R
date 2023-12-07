config_path_dataset <- function(type, acq, region) {
  fs::path(
    "data",
    sprintf("type-%s_acq-%s", type, acq),
    sprintf("region-%s", region)
  )
}
