visualize_mantel <- function(
  patterns_x,
  patterns_y,
  stats,
  name_x,
  name_y,
  show_legend = FALSE
) {
  patterns_flat <- patterns_x |>
    mutate(
      pattern = map(
        pick(last_col())[[1]],
        \(pat) {
          tibble(
            "{name_x}" := unclass(pat),
            "{name_y}" := unclass(patterns_y)
          )
        }
      ),
      .keep = "unused"
    ) |>
    unnest(pattern)
  patterns_flat |>
    ggplot(aes(.data[[name_x]], .data[[name_y]])) +
    geom_hex(
      aes(fill = factor(cca_id), alpha = after_stat(count)),
      show.legend = FALSE
    ) +
    geom_smooth(
      aes(color = factor(cca_id)),
      method = "lm",
      formula = y ~ x,
      show.legend = show_legend
    ) +
    geom_text(
      aes(x = min(patterns_flat[[name_x]]), y = Inf, label = label),
      prepare_corr_plotmath(
        stats,
        "statistic",
        name_p = "italic(p)[Holm]^{Mantel}"
      ),
      hjust = 0,
      vjust = 1,
      parse = TRUE
    ) +
    facet_grid(cols = vars(cca_id), scales = "free") +
    scale_x_continuous(name = name_x) +
    scale_y_continuous(name = name_y) +
    scale_color_components(aesthetics = c("color", "fill")) +
    theme(strip.text = element_blank())
}
