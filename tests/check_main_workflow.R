# Lightweight checks only: no target execution, permutation runs, or report render.
library(tidyverse)
source("R/behav.R")
source("R/utils.R")
source("R_quarto/utils.R")

# Independently check Overall d-prime, including zero and perfect response rates.
events <- tibble(subj_id = rep(1:3, each = 8),
                 old_new = rep(rep(1:2, each = 4), 3),
                 score = c(1, 1, 0, 0, 1, 0, 1, 0,
                           1, 1, 1, 1, 1, 1, 1, 1,
                           0, 0, 0, 0, 0, 0, 0, 0))
expected <- c(0, qnorm(4.5 / 5) - qnorm(0.5 / 5),
              qnorm(0.5 / 5) - qnorm(4.5 / 5))
stopifnot(isTRUE(all.equal(calc_mem_perf(events)$dprime, expected)))

# Statistical wrappers must forward valid arguments and reject misspellings.
data_t <- tibble(cca_id = 1L, time_id = 1L, value = c(1, 3, 2, 5, 4, 6))
stopifnot(calc_stats_t(data_t, value, alternative = "greater")$alternative == "greater")
argument_error <- tryCatch(calc_stats_t(data_t, value, alernative = "greater"),
                           error = identity)
stopifnot(inherits(argument_error, "error"),
          grepl("alernative", conditionMessage(argument_error)))

# Raw significance must not display a cluster that fails Holm correction.
clusters <- tibble(start = c(1, 3, NA), end = c(2, 4, NA),
                   p_perm = c(0.001, 0.03, 1))
annotations <- prepare_cluster_plot_data(clusters)
stopifnot(nrow(annotations) == 1L, annotations$start == 1,
          annotations$p_perm.adj < 0.05)

manifest <- targets::tar_manifest(fields = c(name, command), callr_function = NULL)
stopifnot(!anyDuplicated(manifest$name))
retired <- c("data_combined", "model_med", "fit_med_iss_igs_dprime",
             "fit_med_igs_iss_dprime", "comparison_iss_mem",
             "comparison_iws_mem_charsim_none", "comparison_iss_mem_rt_locked",
             "igs_comp_partial", "iss_comp_partial", "lm_mem_igs_partial",
             "lm_mem_iss_igs_partial", "lm_mem_iss_partial", "lm_mem_igs_iss_partial")
stopifnot(!any(manifest$name %in% retired))

# Parse all authored code without executing Quarto chunks.
for (file in c(list.files("R", "\\.R$", full.names = TRUE),
               list.files("R_quarto", "\\.R$", full.names = TRUE))) {
  invisible(parse(file))
}
purl_path <- tempfile(fileext = ".R")
knitr::purl("index.qmd", output = purl_path, documentation = 0, quiet = TRUE)
report_code <- parse(purl_path)
references <- character()
collect_references <- function(expr) {
  if (!is.call(expr) && !is.expression(expr)) return(invisible(NULL))
  if (is.call(expr) &&
      paste(deparse(expr[[1]]), collapse = "") %in%
      c("targets::tar_read", "targets::tar_load")) {
    references <<- c(references, as.character(expr[[2]]))
  }
  for (child in as.list(expr)) collect_references(child)
}
collect_references(report_code)
missing_targets <- setdiff(references, manifest$name)
stopifnot(length(missing_targets) == 0L,
          !any(grepl("name_perf =", manifest$command, fixed = TRUE)))
cat("PASS:", nrow(manifest), "targets;", length(unique(references)),
    "report target references; local scoring and adjusted cluster annotations.\n")
