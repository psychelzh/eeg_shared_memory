index_time <- function(time_id, onset = 51, sampling_rate = 256) {
  (time_id - onset) / sampling_rate * 1000
}

fit_curve <- function(x, y) {
  nls(
    y ~ eta1 * (1 - exp(theta - eta2 * x)),
    start = list(eta1 = 1, eta2 = 0.01, theta = 0)
  )
}

use_default_name <- function(label = NULL, path = "figures") {
  # https://github.com/rstudio/rstudio/issues/16050
  # unfortunately, opts_current() does not work in interactive mode
  label <- label %||% knitr::opts_current$get("label")
  if (is.null(label) || label == "") {
    stop("No label provided for the figure.")
  }
  fs::path(path, paste0(label, ".png"))
}
