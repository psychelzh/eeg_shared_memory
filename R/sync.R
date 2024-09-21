calc_sync_smc <- function(sync_whole_trials, smc) {
  sync_whole_trials |>
    mutate(
      mantel = map(neu_sync, ~ vegan::mantel(.x, smc)),
      .keep = "unused"
    )
}
