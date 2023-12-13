tar_summary_with_branches <- function(name, data, .by) {
  name <- deparse1(substitute(name))
  list(
    tar_target_raw(
      sprintf("%s_branches", name),
      bquote(
        .(substitute(data)) |>
          lapply(
            summarise,
            n = n(),
            sum_fisher_z = sum(fisher_z),
            .by = .(substitute(.by))
          ) |>
          list_rbind()
      ),
      pattern = bquote(map(.(substitute(data))))
    ),
    tar_target_raw(
      name,
      bquote(
        .(as.symbol(sprintf("%s_branches", name))) |>
          summarise(
            mean_fisher_z = sum(sum_fisher_z) / sum(n),
            .by = .(substitute(.by))
          )
      )
    )
  )
}
