# modified based on https://stackoverflow.com/a/64886536/5996475
visualize_mediation <- function(
  model,
  lab_x,
  lab_y,
  lab_m,
  height = .75,
  width = 2,
  graph_label = "",
  node_text_size = 12,
  edge_text_size = 12,
  color = "black",
  ranksep = .2,
  minlen = 3
) {
  params <- model |>
    parameters::model_parameters(standardize = TRUE) |>
    rstatix::add_significance(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "")
    ) |>
    mutate(
      coef = str_glue("{round(Coefficient, 2)}{p.signif}"),
      style = if_else(p.signif != "", "solid", "dotted")
    ) |>
    as_tibble() |>
    select(Label, coef, style) |>
    pivot_wider(names_from = Label, values_from = c(coef, style))

  # Construct diagram code with Glue
  with(
    params,
    read_file("config/mediation_diagram.dot") |>
      str_glue(.open = "<<", .close = ">>")
  ) |>
    DiagrammeR::grViz()
}
