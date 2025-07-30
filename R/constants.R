# for whole times series analysis, we would remove the first 200 ms baseline
index_onset <- floor(256 * (200 / 1000))
num_subj <- 206L
config_num_subjs <- data.frame(size = seq(20, num_subj, by = 20)) |>
  dplyr::mutate(paired = size <= 100)
config_mem_precise <- tibble::tibble(
  index_name = rlang::syms(c("dprime_rem", "dprime_know"))
)
