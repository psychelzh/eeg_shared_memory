read_eeg_regions <- function(file, channel_regions) {
  data <- qs2::qs_read(file)
  tibble(
    region_id = 1:6,
    data = map(
      region_id,
      \(id) data[channel_regions[[paste0("region", id)]] != 0, , , ]
    )
  )
}

read_charsim <- function(file, mapping, ...) {
  read_tsv(file, show_col_types = FALSE) |>
    pull(similarity) |>
    pracma::squareform() |>
    order_by_trial(mapping)
}

read_rawgray <- function(file, mapping, ...) {
  qs2::qs_read(file) |>
    order_by_trial(mapping)
}

read_alexnet <- function(file, mapping, ..., layer.out = NULL) {
  data <- read_tsv(file, show_col_types = FALSE) |>
    pivot_wider(
      id_cols = c(image_id1, layer),
      names_from = image_id2,
      values_from = similarity
    )
  layer.out <- layer.out %||% unique(data$layer)[[1]]
  data |>
    filter(layer == layer.out) |>
    select(-layer, -image_id1) |>
    as.matrix() |>
    order_by_trial(mapping)
}
