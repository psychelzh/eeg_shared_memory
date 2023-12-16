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

tar_combine_with_meta <- function(name, targets, cols_targets,
                                  prefix = NULL,
                                  fun_pre = NULL,
                                  fun_post = NULL) {
  ischar_name <- tryCatch(
    is.character(name) && length(name) == 1L,
    error = function(e) FALSE
  )
  if (!ischar_name) {
    name <- deparse1(substitute(name))
  }
  if (is.null(prefix)) {
    prefix <- name
  }
  if (is.null(fun_pre)) {
    fun_pre <- \(x) x
  }
  if (is.null(fun_post)) {
    fun_post <- \(x) x
  }
  tarchetypes::tar_combine_raw(
    name,
    targets,
    command = bquote(
      list(!!!.x) |>
        lapply(.(rlang::as_function(fun_pre))) |>
        bind_rows(.id = "id") |>
        # note there is delimiter after prefix should be removed too
        mutate(id = str_remove(id, str_c(.(prefix), "."))) |>
        separate(id, .(cols_targets), convert = TRUE) |>
        .(rlang::as_function(fun_post))()
    )
  )
}
