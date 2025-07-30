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
      coef = str_glue(
        "{round(Coefficient, 2)}{p.signif}"
      ),
      signif = p < 0.05 # Flag for significance (TRUE if p < 0.05, else FALSE)
    ) |>
    as_tibble() |>
    select(Label, coef, signif) |>
    pivot_wider(names_from = Label, values_from = c(coef, signif)) |>
    mutate(
      style_a = ifelse(signif_a, "solid", "dotted"), # x to m
      style_b = ifelse(signif_b, "solid", "dotted"), # m to y
      style_c = ifelse(signif_c, "solid", "dotted") # x to y
    )

  # Construct diagram code with Glue
  diagram_out <- with(
    params,
    str_glue(
      "digraph flowchart {
      fontname = Helvetica
      <<graph_label>>
      graph [ranksep = <<ranksep>>, bgcolor=transparent]

      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle, fixedsize = TRUE, width = <<width>>, height = <<height>>, fontsize = <<node_text_size>>, color = <<color>>]
        mm [label = '<<lab_m>>']
        xx [label = '<<lab_x>>']
        yy [label = '<<lab_y>>']

      # edge definitions with the node IDs
      edge [minlen = <<minlen>>, fontname = Helvetica, fontsize = <<edge_text_size>>, color = <<color>>]
        mm -> yy [label = '<<coef_b>>', style = '<<style_b>>'];
        xx -> mm [label = '<<coef_a>>', style = '<<style_a>>'];
        xx -> yy [label = '<<coef_c>>', style = '<<style_c>>'];

      { rank = same; mm }
      { rank = same; xx; yy }

      }
      ",
      .open = "<<",
      .close = ">>"
    )
  )

  # Generate the diagram with DiagrammeR
  DiagrammeR::grViz(diagram_out)
}
