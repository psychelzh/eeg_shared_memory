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
  ),
  tarchetypes::tar_group_by(
    dat_clean,
    drop_na(dat) |>
      filter(
        n() >= 5,
        .by = c(component, word_category, memory_type, contains("subj_id"))
      ),
    component
  ),
  tar_target(
    fit_sme_clean,
    tibble(
      component = unique(dat_clean$component),
      fit = list(lmerTest::lmer(
        corr_fz ~ word_category * memory_type +
          (1 | subj_id_pair) + (1 | trial_id),
        data = dat_clean
      ))
    ),
    pattern = map(dat_clean)
  ),
  tar_target(
    fit_sme2_clean,
    tibble(
      component = unique(dat_clean$component),
      sapply(
        c("word_category", "memory_type"),
        \(x) {
          lmerTest::lmer(
            str_glue("corr_fz ~ {x} + (1 | subj_id_pair) + (1 | trial_id)"),
            data = dat_clean
          )
        },
        simplify = FALSE
      ) |>
        enframe(name = "term", value = "fit")
    ),
    pattern = map(dat_clean)
  )
)
