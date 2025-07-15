scale_color_components <- function(...) {
  scale_color_manual(
    name = "CCA Comp.",
    values = COLORS_COMPONENTS,
    labels = \(x) paste0("C", x),
    ...
  )
}
