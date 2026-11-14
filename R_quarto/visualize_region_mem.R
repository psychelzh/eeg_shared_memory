visualize_region_mem <- function(
  data_iss,
  data_igs,
  mem_perf,
  stats_iss,
  stats_igs,
  stats_iss_dynamic,
  stats_igs_dynamic,
  clusters_iss_dynamic,
  clusters_igs_dynamic,
  lab_perf = "d'",
  show_uncorrected = TRUE,
  return_panels = FALSE
) {
  region_levels <- paste("Region", 1:6)
  performance_axis <- if (identical(lab_perf, "d'")) {
    expression(italic(d) * "'")
  } else {
    lab_perf
  }
  performance_symbol <- if (identical(lab_perf, "d'")) {
    quote(italic(d) * "'")
  } else {
    lab_perf
  }

  plot_scatter <- function(data, stats, col_stat, lab_stat) {
    data <- data |>
      left_join(mem_perf, by = "subj_id") |>
      filter(is.finite(dprime)) |>
      mutate(region = factor(region_id, levels = 1:6, labels = region_levels))
    stats <- stats |>
      prepare_corr_plotmath() |>
      mutate(
        region = factor(region_id, levels = 1:6, labels = region_levels),
        x = min(data[[col_stat]]),
        y = Inf
      )
    data |>
      ggplot(aes(.data[[col_stat]], dprime)) +
      geom_point(color = "grey40", alpha = 0.7) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        color = "black",
        linewidth = 0.7
      ) +
      geom_text(
        data = stats,
        aes(x, y, label = label),
        hjust = 0,
        vjust = 1,
        parse = TRUE,
        inherit.aes = FALSE
      ) +
      facet_wrap(vars(region), ncol = 3, scales = "free_x") +
      scale_x_continuous(name = lab_stat) +
      scale_y_continuous(name = performance_axis)
  }

  plot_dynamic <- function(stats, clusters, lab_stat) {
    stats <- stats |>
      mutate(
        region = factor(region_id, levels = 1:6, labels = region_levels),
        time = index_time(time_id)
      )
    limits <- range(c(stats$conf.low, stats$conf.high))
    span <- max(diff(limits), .Machine$double.eps)
    pointwise <- stats |>
      filter(is.finite(p.value), p.value < 0.05)
    exploratory_clusters <- clusters |>
      ungroup() |>
      mutate(p_holm = p.adjust(p_perm, method = "holm")) |>
      filter(!is.na(start), !is.na(end), p_perm < 0.05, p_holm >= 0.05) |>
      mutate(
        region = factor(region_id, levels = 1:6, labels = region_levels),
        label = sprintf("italic(p)[cl] == %.3f", p_perm)
      )
    corrected_clusters <- clusters |>
      prepare_cluster_plot_data() |>
      mutate(region = factor(region_id, levels = 1:6, labels = region_levels))
    plot_limits <- if (show_uncorrected) {
      limits + c(-0.10, 0.12) * span
    } else {
      limits
    }
    stats |>
      ggplot(aes(time, estimate)) +
      {
        if (show_uncorrected && nrow(exploratory_clusters) > 0) {
          geom_rect(
            data = exploratory_clusters,
            aes(
              xmin = index_time(start),
              xmax = index_time(end),
              ymin = limits[1],
              ymax = limits[2]
            ),
            inherit.aes = FALSE,
            fill = "#E69F00",
            colour = "#A65B00",
            alpha = 0.12,
            linetype = "dashed",
            linewidth = 0.4
          )
        }
      } +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
      geom_line(linewidth = 0.8) +
      {
        if (nrow(corrected_clusters) > 0) {
          list(
            geom_rect(
              data = corrected_clusters,
              aes(
                xmin = index_time(start),
                xmax = index_time(end),
                ymin = limits[1],
                ymax = limits[2]
              ),
              inherit.aes = FALSE,
              alpha = 0.1
            ),
            geom_text(
              data = corrected_clusters,
              aes(
                x = index_time((start + end) / 2),
                y = limits[2],
                label = p_perm.adj.signif
              ),
              size = SIZE_LABEL,
              inherit.aes = FALSE
            )
          )
        }
      } +
      {
        if (show_uncorrected && nrow(pointwise) > 0) {
          geom_segment(
            data = pointwise,
            aes(
              x = time,
              xend = time,
              y = limits[1] - 0.07 * span,
              yend = limits[1] - 0.025 * span
            ),
            inherit.aes = FALSE,
            colour = "#0072B2",
            linewidth = 0.7
          )
        }
      } +
      {
        if (show_uncorrected && nrow(exploratory_clusters) > 0) {
          geom_text(
            data = exploratory_clusters,
            aes(
              x = index_time((start + end) / 2),
              y = limits[2] + 0.065 * span,
              label = label
            ),
            inherit.aes = FALSE,
            parse = TRUE,
            colour = "#A65B00",
            size = 2.8
          )
        }
      } +
      facet_wrap(vars(region), ncol = 3) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey") +
      scale_x_continuous(name = "Encoding Time (ms)") +
      scale_y_continuous(name = lab_stat, limits = plot_limits)
  }

  p_iss_scatter <- plot_scatter(data_iss, stats_iss, "iss", "ISS")
  p_igs_scatter <- plot_scatter(data_igs, stats_igs, "igs", "IGS")
  p_iss_dynamic <- plot_dynamic(
    stats_iss_dynamic,
    clusters_iss_dynamic,
    bquote(italic(r)[ISS - .(performance_symbol)])
  )
  p_igs_dynamic <- plot_dynamic(
    stats_igs_dynamic,
    clusters_igs_dynamic,
    bquote(italic(r)[IGS - .(performance_symbol)])
  )

  has_corrected_cluster <- any(vapply(
    list(clusters_iss_dynamic, clusters_igs_dynamic),
    function(x) {
      any(!is.na(x$start) & !is.na(x$end) & p.adjust(x$p_perm, "holm") < 0.05)
    },
    logical(1)
  ))
  holm_note <- if (has_corrected_cluster) {
    "Grey bands/stars: Holm-significant positive clusters; negative clusters were not tested."
  } else {
    "No positive cluster survives Holm; negative clusters were not tested."
  }
  caption <- if (show_uncorrected) {
    paste(
      "Blue ticks: pointwise significance (two-sided, uncorrected; threshold .05).",
      "Orange dashed bands: positive, one-sided cluster significance before Holm only (threshold .05).",
      holm_note,
      sep = "\n"
    )
  }
  panels <- list(p_iss_scatter, p_igs_scatter, p_iss_dynamic, p_igs_dynamic)
  attr(panels, "significance_caption") <- caption
  if (return_panels) {
    return(panels)
  }
  (panels[[1]] | panels[[2]]) /
    (panels[[3]] | panels[[4]]) +
    plot_annotation(
      tag_levels = "A",
      caption = caption,
      theme = theme(plot.caption = element_text(hjust = 0, size = 9.5))
    )
}
