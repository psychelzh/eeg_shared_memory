library(targets)

tar_option_set(
  packages = c("tidyverse"),
  format = "qs",
  controller = crew::crew_controller_local(
    name = "local",
    workers = 4
  )
)

list(
  tar_target(
    file_isc,
    "data/CorCA/isc_across_time_each_trial_first3components.parquet",
    format = "file"
  ),
  tar_target(
    file_mem,
    "data/CorCA/cross-subject_task-wordencoding_labels_pairwise.parquet",
    format = "file"
  ),
  tarchetypes::tar_group_by(
    dat,
    arrow::read_parquet(file_mem) |>
      filter(if_all(contains("memory_type"), \(x) x != 0)) |>
      filter(i1_memory_type == i2_memory_type) |>
      mutate(memory_type = factor(i1_memory_type)) |>
      left_join(
        arrow::read_parquet(file_isc),
        by = join_by(i1_subj_id, i2_subj_id, trial_id)
      ) |>
      mutate(subj_id_pair = paste(i1_subj_id, i2_subj_id, sep = "_")),
    component
  ),
  tar_target(
    fit_sme,
    tibble(
      component = unique(dat$component),
      fit = list(lmerTest::lmer(
        corr_fz ~ word_category * memory_type +
          (word_category * memory_type | subj_id_pair) + (1 | trial_id),
        data = dat
      ))
    ),
    pattern = map(dat)
  )
)
