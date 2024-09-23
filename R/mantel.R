calc_mantel <- function(data, ydis, col_xdis = last_col()) {
  data |>
    mutate(
      mantel = map(
        pick({{ col_xdis }})[[1]],
        \(xdis) vegan::mantel(xdis, ydis)
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
