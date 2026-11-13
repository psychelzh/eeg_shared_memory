# Run after tar_make() to export statistical summaries from the completed cache.
# This script reads targets; it does not rebuild or invalidate the pipeline.
library(tidyverse)
library(targets)
source("R_quarto/utils.R")

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) {
  args[[1L]]
} else {
  file.path("docs", paste0("analysis-update-", Sys.Date()))
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
manifest <- tar_manifest(fields = name, callr_function = NULL)

write_summary <- function(x, filename) {
  readr::write_csv(x, file.path(output, filename), na = "")
}

flatten <- function(x) {
  if ("cluster_mass_perm" %in% names(x)) {
    x$n_permutations <- lengths(x$cluster_mass_perm)
  }
  x[, !vapply(x, is.list, logical(1)), drop = FALSE]
}

whole_names <- manifest$name[str_detect(manifest$name, "^stats_.*whole")]
whole <- map(whole_names, function(name) {
  x <- flatten(tar_read_raw(name))
  if ("p.value" %in% names(x)) {
    x$p_holm <- p.adjust(x$p.value, "holm")
  }
  if ("p_perm" %in% names(x)) {
    x$p_perm_holm <- p.adjust(x$p_perm, "holm")
  }
  add_column(x, target = name, .before = 1)
}) |>
  list_rbind()
write_summary(whole, "whole_time_statistics.csv")

cluster_names <- manifest$name[str_detect(manifest$name, "^clusters_stats_")]
clusters <- map(cluster_names, function(name) {
  x <- flatten(tar_read_raw(name)) |>
    mutate(
      p_holm = p.adjust(p_perm, "holm"),
      significant_holm = !is.na(start) & !is.na(end) & p_holm < 0.05,
      start_ms = index_time(start),
      end_ms = index_time(end)
    )
  add_column(x, target = name, .before = 1)
}) |>
  list_rbind()
write_summary(clusters, "temporal_clusters.csv")
write_summary(
  filter(clusters, significant_holm),
  "significant_temporal_clusters.csv"
)

memory <- tar_read(mem_perf) |>
  left_join(tar_read(mem_perf_precise), by = "subj_id")
behavior <- map(setdiff(names(memory), "subj_id"), function(name) {
  psych::describe(memory[[name]]) |>
    as_tibble() |>
    mutate(measure = name, .before = 1)
}) |>
  list_rbind()
write_summary(behavior, "memory_descriptives.csv")
write_summary(
  as_tibble(cor(select(memory, -subj_id)), rownames = "measure"),
  "memory_correlations.csv"
)
extra_behavior <- tibble(
  measure = c(
    "SMC",
    "memorability_content_r",
    "memorability_content_r_squared"
  ),
  mean = c(
    mean(tar_read(smc)),
    mean(tar_read(memorability_content)$r),
    mean(tar_read(memorability_content)$r^2)
  ),
  sd = c(
    sd(tar_read(smc)),
    sd(tar_read(memorability_content)$r),
    sd(tar_read(memorability_content)$r^2)
  )
)
write_summary(extra_behavior, "additional_behavior_descriptives.csv")

write_summary(
  map(c("iss", "igs", "gss", "isps"), function(measure) {
    x <- tar_read_raw(paste0("stats_", measure, "_dynamic"))
    column <- if (measure == "isps") {
      "isps_mean"
    } else if (measure == "gss") {
      "statistic"
    } else {
      "estimate"
    }
    x |>
      slice_max(.data[[column]], n = 1, with_ties = FALSE, by = cca_id) |>
      transmute(
        measure,
        cca_id,
        time_id,
        peak_ms = index_time(time_id),
        value = .data[[column]]
      )
  }) |>
    list_rbind(),
  "temporal_peaks.csv"
)

comparison_names <- intersect(
  c("iss_comparison", "igs_comparison"),
  manifest$name
)
write_summary(
  map(comparison_names, function(name) {
    add_column(flatten(tar_read_raw(name)), target = name, .before = 1)
  }) |>
    list_rbind(),
  "component_comparisons.csv"
)

joint <- tar_read(data_combined_summary_sum)
fit <- joint |>
  transmute(across(c(dprime, c_iss, c_igs), ~ as.numeric(scale(.x)))) |>
  lm(dprime ~ c_iss + c_igs, data = _)
write_summary(broom::glance(fit), "joint_model.csv")
write_summary(broom::tidy(fit, conf.int = TRUE), "joint_model_coefficients.csv")
joint_correlations <- tribble(
  ~x,
  ~y,
  "c_iss",
  "c_igs",
  "c_iss",
  "dprime",
  "c_igs",
  "dprime"
) |>
  mutate(
    test = map2(x, y, ~ broom::tidy(cor.test(joint[[.x]], joint[[.y]])))
  ) |>
  unnest(test)
write_summary(joint_correlations, "joint_score_correlations.csv")
write_summary(tar_read(ctb_model_iss_igs_dprime), "commonality.csv")
write_summary(
  tar_read(patterns_group_stability) |>
    summarise(mean_r = mean(r), sd_r = sd(r), .by = c(cca_id, size)),
  "group_stability.csv"
)

# Basic export integrity checks; plot and manuscript corrections use the same scope.
stopifnot(
  nrow(memory) == 206L,
  all(clusters$n_permutations == 1000L),
  all(clusters$p_holm >= 0 & clusters$p_holm <= 1),
  all(whole$p_holm >= 0 & whole$p_holm <= 1, na.rm = TRUE)
)
cat(
  "Exported",
  length(whole_names),
  "whole-time targets and",
  length(cluster_names),
  "cluster targets to",
  output,
  "\n"
)
