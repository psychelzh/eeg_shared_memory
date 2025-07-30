index_time <- function(time_id, onset = 51, sampling_rate = 256) {
  (time_id - onset) / sampling_rate * 1000
}

fit_curve <- function(x, y) {
  nls(
    y ~ eta1 * (1 - exp(theta - eta2 * x)),
    start = list(eta1 = 1, eta2 = 0.01, theta = 0)
  )
}

ggsave_default <- purrr::partial(
  ggplot2::ggsave,
  filename = use_default_name(),
  dpi = 600,
  bg = "transparent"
)

use_default_name <- function(label = NULL, path = "figures") {
  # https://github.com/rstudio/rstudio/issues/16050
  # unfortunately, opts_current() does not work in interactive mode
  label <- label %||% knitr::opts_current$get("label")
  if (is.null(label) || label == "") {
    stop("No label provided for the figure.")
  }
  fs::path(path, paste0(label, ".png"))
}

prepare_corr_plotmath <- function(
  stats,
  col_r = "estimate",
  col_p = "p.value",
  name_r = "italic(r)",
  name_p = "italic(p)[Holm]"
) {
  stats |>
    rstatix::adjust_pvalue(col_p, "p_adj") |>
    rstatix::add_significance(
      "p_adj",
      "p_adj_sig",
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "")
    ) |>
    mutate(
      label = format_r_plotmath(
        .data[[col_r]],
        p_adj,
        p.sig = p_adj_sig,
        name_r = name_r,
        name_p = name_p
      )
    )
}

format_r_plotmath <- function(
  r,
  p,
  p.sig = "",
  name_r = "italic(r)",
  name_p = "italic(p)[Holm]"
) {
  paste0(
    str_glue("{name_r}*' = '*{round(r, 2)}"),
    if (is.null(name_p)) {
      str_glue("^'{p.sig}'")
    } else {
      paste0(
        "*', '*",
        if_else(
          p < 0.001,
          str_glue("{name_p} < 0.001^'{p.sig}'"),
          str_glue("{name_p}*' = '*{round(p, 3)}^'{p.sig}'")
        )
      )
    }
  )
}
