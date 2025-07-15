visualize_dynamic <- function(
  stats,
  clusters_stats = NULL,
  col_stat = "estimate",
  lab_stat = "Estimate",
  col_cis = c("conf.low", "conf.high"),
  limits = NULL,
  show_legend = FALSE
) {
  if (!is.null(clusters_stats)) {
    clusters_stats <- clusters_stats |>
      mutate(cca_id = factor(cca_id)) |>
      rstatix::adjust_pvalue("p_perm") |>
      rstatix::add_significance(
        "p_perm.adj",
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "")
      ) |>
      filter(p_perm < 0.05)
  }
  show_cis <- !is.null(col_cis) && all(has_name(stats, col_cis))
  limits_rect <- if (show_cis) {
    range(c(stats[[col_cis[1]]], stats[[col_cis[2]]]))
  } else {
    range(stats[[col_stat]])
  }
  stats |>
    mutate(
      cca_id = factor(cca_id),
      time = index_time(time_id)
    ) |>
    ggplot(aes(time, .data[[col_stat]])) +
    geom_line(
      aes(color = cca_id),
      linewidth = 1,
      show.legend = show_legend
    ) +
    {
      # use ribbon to show confidence intervals if available
      if (show_cis) {
        geom_ribbon(
          aes(
            fill = cca_id,
            ymin = .data[[col_cis[1]]],
            ymax = .data[[col_cis[2]]]
          ),
          alpha = 0.2,
          show.legend = show_legend
        )
      }
    } +
    {
      # annotate significant clusters using rect and text
      if (!is.null(clusters_stats)) {
        list(
          geom_rect(
            data = clusters_stats,
            mapping = aes(
              xmin = index_time(start),
              xmax = index_time(end),
              ymin = limits_rect[1],
              ymax = limits_rect[2]
            ),
            inherit.aes = FALSE,
            alpha = 0.1
          ),
          geom_text(
            data = clusters_stats,
            mapping = aes(
              x = index_time((start + end) / 2),
              y = limits_rect[2],
              label = p_perm.adj.signif
            ),
            size = SIZE_LABEL,
            inherit.aes = FALSE
          )
        )
      }
    } +
    facet_grid(cols = vars(cca_id)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey") +
    scale_x_continuous(name = "Encoding Time (ms)") +
    scale_y_continuous(name = lab_stat, limits = limits) +
    scale_color_components(aesthetics = c("color", "fill")) +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank()
    )
}
