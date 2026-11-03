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
  lab_perf = "d'"
) {
  region_levels <- paste("Region", 1:6)

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
      scale_y_continuous(name = lab_perf)
  }

  plot_dynamic <- function(stats, clusters, lab_stat) {
    stats <- stats |>
      mutate(
        region = factor(region_id, levels = 1:6, labels = region_levels),
        time = index_time(time_id)
      )
    limits <- range(c(stats$conf.low, stats$conf.high))
    clusters <- clusters |>
      prepare_cluster_plot_data() |>
      mutate(region = factor(region_id, levels = 1:6, labels = region_levels))
    stats |>
      ggplot(aes(time, estimate)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
      geom_line(linewidth = 0.8) +
      {
        if (nrow(clusters) > 0) {
          list(
            geom_rect(
              data = clusters,
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
              data = clusters,
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
      facet_wrap(vars(region), ncol = 3) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey") +
      geom_vline(xintercept = 0, linetype = "dotted", color = "grey") +
      scale_x_continuous(name = "Encoding Time (ms)") +
      scale_y_continuous(name = lab_stat, limits = limits)
  }

  p_iss_scatter <- plot_scatter(data_iss, stats_iss, "iss", "ISS")
  p_igs_scatter <- plot_scatter(data_igs, stats_igs, "igs", "IGS")
  p_iss_dynamic <- plot_dynamic(
    stats_iss_dynamic,
    clusters_iss_dynamic,
    bquote(italic(r)[ISS - .(lab_perf)])
  )
  p_igs_dynamic <- plot_dynamic(
    stats_igs_dynamic,
    clusters_igs_dynamic,
    bquote(italic(r)[IGS - .(lab_perf)])
  )

  (p_iss_scatter | p_igs_scatter) /
    (p_iss_dynamic | p_igs_dynamic) +
    plot_annotation(tag_levels = "A")
}
