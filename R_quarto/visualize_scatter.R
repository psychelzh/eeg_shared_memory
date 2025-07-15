visualize_scatter <- function(
  data,
  mem_perf,
  lab_stat,
  col_stat,
  show_legend = FALSE
) {
  data_joind <- data |>
    left_join(mem_perf, by = "subj_id") |>
    mutate(cca_id = factor(cca_id))
  stats <- data_joind |>
    reframe(
      cor.test(.data[[col_stat]], .data$dprime) |>
        broom::tidy(),
      .by = cca_id
    ) |>
    prepare_corr_plotmath()
  data_joind |>
    ggplot(aes(.data[[col_stat]], dprime)) +
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
    scale_y_continuous(name = "d'") +
    scale_color_components() +
    theme(strip.text = element_blank())
}
