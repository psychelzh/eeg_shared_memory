visualize_scatter <- function(
  data,
  mem_perf,
  col_stat = NULL,
  col_perf = NULL,
  lab_stat = NULL,
  lab_perf = NULL,
  show_legend = FALSE
) {
  # use the last column as default if col_stat is not provided
  col_stat <- col_stat %||% last(names(data))
  col_perf <- col_perf %||% last(names(mem_perf))
  lab_stat <- lab_stat %||% str_to_upper(col_stat)
  lab_perf <- lab_perf %||% "d'"
  data_joind <- data |>
    left_join(mem_perf, by = "subj_id") |>
    mutate(cca_id = factor(cca_id))
  stats <- data_joind |>
    reframe(
      cor.test(.data[[col_stat]], .data[[col_perf]]) |>
        broom::tidy(),
      .by = cca_id
    ) |>
    prepare_corr_plotmath()
  data_joind |>
    ggplot(aes(.data[[col_stat]], .data[[col_perf]])) +
    geom_point(aes(color = cca_id), show.legend = show_legend) +
    geom_smooth(
      aes(color = cca_id),
      method = "lm",
      formula = y ~ x,
      show.legend = show_legend
    ) +
    geom_text(
      aes(x = min(data_joind[[col_stat]]), y = Inf, label = label),
      stats,
      hjust = 0,
      vjust = 1,
      parse = TRUE
    ) +
    facet_grid(cols = vars(cca_id), scales = "free") +
    scale_x_continuous(name = lab_stat) +
    scale_y_continuous(name = lab_perf) +
    scale_color_components() +
    theme(strip.text = element_blank())
}
