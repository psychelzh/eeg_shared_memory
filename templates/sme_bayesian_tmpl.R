library(tidyverse)
library(brms)
file_isc <- "data/CorCA/isc_across_time_each_trial_first3components.csv"
file_mem <- "data/CorCA/cross-subject_task-wordencoding_labels_pairwise.csv"
read_csv(file_mem, show_col_types = FALSE) |>
  filter(if_all(contains("memory_type"), \(x) x != 0)) |>
  filter(i1_memory_type == i2_memory_type) |>
  mutate(memory_type = factor(i1_memory_type)) |>
  left_join(
    read_csv(file_isc, show_col_types = FALSE),
    by = join_by(i1_subj_id, i2_subj_id, trial_id)
  ) |>
  mutate(subj_id_pair = paste(i1_subj_id, i2_subj_id, sep = "_")) |>
  # fmt: skip
  filter(component == {component}) |>
  brm(
    formula = corr_fz ~
      0 +
        word_category:memory_type +
        (0 + word_category:memory_type | subj_id_pair) +
        (1 | trial_id),
    data = _,
    chains = 4,
    cores = 4,
    iter = 2000,
    warmup = 1000,
    file = file.path(
      "results/",
      paste0("comp-{component}_fitSME")
    )
  )
