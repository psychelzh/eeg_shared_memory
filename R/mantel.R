calc_mantel <- function(data, ydis, col_xdis = last_col()) {
  data |>
    mutate(
      mantel = map(
        pick({{ col_xdis }})[[1]],
        \(xdis) vegan::mantel(xdis, ydis, permutations = 9999)
      ),
      .keep = "unused"
    )
}

calc_mantel_partial <- function(data, ydis, zdis, col_xdis = last_col()) {
  data |>
    mutate(
      mantel = map(
        pick({{ col_xdis }})[[1]],
        \(xdis) vegan::mantel.partial(xdis, ydis, zdis, permutations = 9999)
      ),
      .keep = "unused"
    )
}

calc_mantel2 <- function(
  data_xdis,
  data_ydis
) {
  data_xdis |>
    cross_join(data_ydis) |>
    mutate(
      mantel = map2(
        pattern.x,
        pattern.y,
        \(xdis, ydis) vegan::mantel(xdis, ydis, permutations = 9999)
      ),
      .keep = "unused"
    )
}

extract_stats_mantel <- function(data_mantel) {
  data_mantel |>
    mutate(
      map(mantel, tidy_mantel) |>
        list_rbind(),
      .keep = "unused"
    )
}

tidy_mantel <- function(mantel) {
  tibble(
    statistic = mantel$statistic,
    p.value = mantel$signif,
    method = mantel$method
  )
}
