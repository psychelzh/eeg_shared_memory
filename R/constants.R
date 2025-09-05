# for whole times series analysis, we would remove the first 200 ms baseline
index_onset <- floor(256 * (200 / 1000))
num_subj <- 206L
config_num_subjs <- data.frame(size = seq(20, num_subj, by = 20)) |>
  dplyr::mutate(paired = size <= 100)
config_mem_precise <- tibble::tibble(
  index_name = rlang::syms(c("dprime_rem", "dprime_know"))
)
word_shape_methods <- dplyr::bind_rows(
  tibble::tibble(model = c("charsim", "rawgray"), layer = "none"),
  tibble::tibble(
    model = "alexnet",
    layer = c(
      "conv2d_1_1",
      "conv2d_2_4",
      "conv2d_3_7",
      "conv2d_4_9",
      "conv2d_5_11",
      "linear_1_17",
      "linear_2_20",
      "linear_3_22"
    )
  )
) |>
  dplyr::mutate(
    file = rlang::syms(paste0("file_", model)),
    read = rlang::syms(paste0("read_", model)),
    name_pattern = rlang::syms(
      paste("pattern_shapes", model, layer, sep = "_")
    )
  )
