visualize_comp_preds <- function(...) {
  list_preds <- list(...)
  preds <- list_preds |>
    map("pred") |>
    list_rbind(names_to = "model") |>
    mutate(model = factor(model, names(list_preds)))
  model_eval <- tibble(
    x = min(preds$obs),
    y = max(preds$pred) * (1 + c(0.1, 0.02)),
    model = factor(names(list_preds)),
    rsquared = map_dbl(
      list_preds,
      \(x) caret::getTrainPerf(x)$TrainRsquared
    ),
    label = str_c(
      "*R*<sup>2</sup><sub>",
      model,
      "</sub> = ",
      signif(rsquared, 2)
    )
  )
  preds |>
    ggplot(aes(obs, pred, color = model)) +
    geom_point(shape = 16) +
    geom_smooth(method = "lm", formula = y ~ x) +
    ggtext::geom_richtext(
      aes(x, y, color = model, label = label),
      model_eval,
      size = 3,
      fill = NA,
      label.color = NA, # remove background and outline
      label.padding = grid::unit(rep(0, 4), "pt"), # remove padding
      hjust = 0,
      vjust = 0.5, # bottom-left corner
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    # facet_grid(cols = vars(model)) +
    scale_x_continuous(name = "Observed") +
    scale_y_continuous(name = "Predicted") +
    scale_color_grey(start = 0.1, end = 0.6, name = "Model", guide = "none") +
    theme(strip.background = element_blank())
}
