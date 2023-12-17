tar_combine_with_meta <- function(name, targets, cols_targets, ...,
                                  prefix = NULL,
                                  fun_pre = NULL,
                                  fun_post = NULL) {
  rlang::check_dots_used()
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
    ),
    ...
  )
}
