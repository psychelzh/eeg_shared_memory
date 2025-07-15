visualize_mantel_dist <- function(data, stats, label, show_legend = FALSE) {
  data |>
    mutate(
      cca_id = factor(cca_id),
      null = map(mantel, "perm"),
      .keep = "unused"
    ) |>
    unchop(null) |>
    ggplot(aes(null)) +
    geom_histogram(aes(fill = cca_id), show.legend = show_legend) +
    geomtextpath::geom_textvline(
      aes(xintercept = statistic, label = label),
      stats |>
        mutate(cca_id = factor(cca_id)) |>
        prepare_corr_plotmath(
          "statistic",
          name_r = "italic(r)[Obs]",
          name_p = NULL
        ),
      parse = TRUE,
      vjust = -0.1
    ) +
    facet_grid(cols = vars(cca_id)) +
    scale_x_continuous(name = label) +
    scale_y_continuous(name = "Count", expand = expansion(mult = c(0, 0.05))) +
    scale_color_components(aesthetics = "fill") +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank()
    )
}
