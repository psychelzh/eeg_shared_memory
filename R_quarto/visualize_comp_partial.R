visualize_comp_partial <- function(
  comp_partial,
  ...,
  mult = 1,
  comp_labels = waiver()
) {
  stats <- comp_partial$stats |>
    mutate(
      ymax = estimate + mult * std.error,
      ymin = estimate - mult * std.error
    )
  htests <- comp_partial$htest |>
    separate_wider_delim(
      contrast,
      " - ",
      names = c("start", "end")
    ) |>
    rstatix::add_significance(
      "p.value",
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "")
    ) |>
    mutate(y_position = max(stats$ymax) * 1.05) |>
    filter(p.value < 0.05)
  stats |>
    ggplot(aes(type, estimate, ymax = ymax, ymin = ymin, color = cca_id)) +
    geom_point(size = SIZE_LABEL) +
    geom_errorbar(width = 0.1, linewidth = 1) +
    geom_line(aes(group = cca_id), linewidth = 1) +
    ggsignif::geom_signif(
      aes(
        xmin = start,
        xmax = end,
        annotations = p.value.signif,
        y_position = y_position
      ),
      htests,
      size = 0.8,
      textsize = SIZE_LABEL,
      vjust = 0.5,
      inherit.aes = FALSE,
      manual = TRUE
    ) +
    facet_grid(cols = vars(cca_id), scales = "free") +
    scale_x_discrete(name = NULL, labels = comp_labels) +
    scale_y_continuous(name = expression(italic(r))) +
    scale_color_components(guide = "none") +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank()
    )
}
