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
    theme(axis.line = element_line(linewidth = 1), strip.text = element_blank())
}

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
    theme(
      strip.text = element_blank(),
      axis.line = element_line(linewidth = 1)
    )
}

visualize_mantel_dist <- function(data, stats, label, show_legend = FALSE) {
  data |>
    mutate(
      cca_id = factor(cca_id),
      null = map(mantel, "perm"),
      .keep = "unused"
    ) |>
    unchop(null) |>
    ggplot(aes(null)) +
    geom_histogram(aes(fill = cca_id), show.legend = show_legend) +
    geomtextpath::geom_textvline(
      aes(xintercept = statistic, label = label),
      stats |>
        mutate(cca_id = factor(cca_id)) |>
        prepare_corr_plotmath(
          "statistic",
          name_r = "italic(r)[Obs]",
          name_p = NULL
        ),
      parse = TRUE,
      vjust = -0.1
    ) +
    facet_grid(cols = vars(cca_id)) +
    scale_x_continuous(name = label) +
    scale_y_continuous(name = "Count", expand = expansion(mult = c(0, 0.05))) +
    scale_color_components(aesthetics = "fill") +
    theme(
      axis.line = element_line(linewidth = 1),
      strip.text = element_blank(),
      strip.background = element_blank()
    )
}

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
    # TODO: Convert these as functions
    {
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
            size = size_label,
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
      strip.background = element_blank(),
      axis.line = element_line(linewidth = 1)
    )
}

scale_color_components <- function(...) {
  scale_color_manual(
    name = "CCA Comp.",
    values = colors_components,
    labels = \(x) paste0("C", x),
    ...
  )
}

# internal functions
prepare_corr_plotmath <- function(
  stats,
  col_r = "estimate",
  col_p = "p.value",
  name_r = "italic(r)",
  name_p = "italic(p)[Holm]"
) {
  stats |>
    rstatix::adjust_pvalue(col_p, "p_adj") |>
    rstatix::add_significance(
      "p_adj",
      "p_adj_sig",
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "")
    ) |>
    mutate(
      label = format_r_plotmath(
        .data[[col_r]],
        p_adj,
        p.sig = p_adj_sig,
        name_r = name_r,
        name_p = name_p
      )
    )
}

format_r_plotmath <- function(
  r,
  p,
  p.sig = "",
  name_r = "italic(r)",
  name_p = "italic(p)[Holm]"
) {
  paste0(
    str_glue("{name_r}*' = '*{round(r, 2)}"),
    if (is.null(name_p)) {
      str_glue("^'{p.sig}'")
    } else {
      paste0(
        "*', '*",
        if_else(
          p < 0.001,
          str_glue("{name_p} < 0.001^'{p.sig}'"),
          str_glue("{name_p}*' = '*{round(p, 3)}^'{p.sig}'")
        )
      )
    }
  )
}
