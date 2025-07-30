visualize_data_distr <- function(
  data,
  stats,
  comparison,
  ...,
  col_data = last_col(),
  col_label = NULL,
  show_legend = FALSE,
  signif_base = 0.03,
  signif_step = 0.12
) {
  rlang::check_dots_empty()
  col_data_name <- names(select(data, {{ col_data }}))
  col_label <- col_label %||% str_to_upper(col_data_name)
  max_y <- max(data[[col_data_name]], na.rm = TRUE)
  stats <- stats |>
    rstatix::adjust_pvalue() |>
    rstatix::add_significance(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "")
    ) |>
    add_column(y = max_y)
  comparison <- comparison |>
    filter(adj.p.value < 0.05) |>
    mutate(
      across(c(start, end), \(x) factor(x, levels = 1:3)),
      y_position = max_y * (1 + signif_base + signif_step * seq_len(n()))
    ) |>
    rstatix::add_significance(
      "adj.p.value",
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "")
    )
  data |>
    mutate(cca_id = factor(cca_id)) |>
    ggplot(aes(cca_id, .data[[col_data_name]])) +
    ggdist::stat_dotsinterval(
      aes(slab_color = cca_id, slab_fill = cca_id),
      slab_alpha = 0.4,
      side = "both"
    ) +
    geom_text(
      aes(cca_id, y, label = p.value.adj.signif),
      stats,
      size = SIZE_LABEL,
      vjust = 0,
      inherit.aes = FALSE
    ) +
    ggsignif::geom_signif(
      aes(
        xmin = start,
        xmax = end,
        annotations = adj.p.value.signif,
        y_position = y_position
      ),
      comparison,
      size = 0.8,
      textsize = SIZE_LABEL,
      vjust = 0.5,
      inherit.aes = FALSE,
      manual = TRUE
    ) +
    scale_x_discrete(name = NULL, labels = \(x) paste0("C", x)) +
    scale_y_continuous(name = col_label) +
    scale_color_components(
      aesthetics = c("slab_color", "slab_fill"),
      guide = if (show_legend) "legend" else "none"
    )
}
