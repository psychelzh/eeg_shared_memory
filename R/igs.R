calc_igs <- function(patterns, patterns_group) {
  patterns |>
    mutate(
      igs = map2_dbl(
        pattern, cca_id,
        \(pat, cca) {
          with(
            patterns_group,
            cor(
              pattern[[which(cca_id == cca)]],
              pat,
              use = "pairwise"
            )
          )
        }
      )
    ) |>
    select(!pattern)
}
