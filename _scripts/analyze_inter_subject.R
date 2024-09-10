library(targets)

tar_option_set(
  packages = c("tidyverse"),
  controller = crew::crew_controller_local(
    name = "local",
    workers = 8
  )
)

list(
  tar_target(file_cca_y, "data/CorCAExtra/cca_y_subjs206_flat.parquet"),
  tar_target(
    subj_id_loop,
    arrow::open_dataset(file_cca_y) |>
      distinct(subj_id) |>
      pull(subj_id, as_vector = TRUE)
  ),
  tar_target(
    patterns_cca,
    arrow::open_dataset(file_cca_y) |>
      filter(subj_id == subj_id_loop) |>
      collect() |>
      pivot_wider(names_from = trial_id, values_from = y) |>
      reframe(
        pick(!time_id) |>
          slider::slide(
            \(x) as.dist(cor(x, use = "pairwise")),
            .before = 25,
            .after = 25,
            .step = 5,
            .complete = TRUE
          ) |>
          enframe(name = "time_id", value = "pattern_cca") |>
          filter(!map_lgl(pattern_cca, is.null)),
        .by = c(subj_id, cca_id)
      ),
    pattern = map(subj_id_loop),
    iteration = "list"
  )
)
