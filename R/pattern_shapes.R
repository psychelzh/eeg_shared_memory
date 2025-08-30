read_charsim <- function(file, mapping) {
  read_tsv(file, show_col_types = FALSE) |>
    pull(similarity) |>
    pracma::squareform() |>
    order_by_trial(mapping)
}

read_alexnet <- function(file, mapping) {
  read_tsv(file, show_col_types = FALSE) |>
    pivot_wider(
      id_cols = c(image_id1, layer),
      names_from = image_id2,
      values_from = similarity
    ) |>
    nest(.by = layer) |>
    mutate(
      pattern = map(data, \(x) {
        x |>
          select(-image_id1) |>
          as.matrix() |>
          order_by_trial(mapping)
      }),
      .keep = "unused"
    ) |>
    add_column(model = "alexnet", .before = 1)
}

order_by_trial <- function(pattern, mapping) {
  order <- with(mapping, word_id[trial_id > 0])
  as.dist(pattern[order, order])
}
