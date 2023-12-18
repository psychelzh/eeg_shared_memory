calc_dist <- function(dat, fun) {
  dat_mat <- as.matrix(dat)
  size <- nrow(dat)
  mat <- matrix(0, nrow = size, ncol = size)
  for (i_col in seq_len(size - 1)) {
    for (i_row in seq(i_col + 1, size)) {
      mat[i_row, i_col] <-
        fun(dat_mat[i_row, ], dat_mat[i_col, ])
    }
  }
  rownames(mat) <- rownames(dat)
  as.dist(mat)
}

calc_dist_mem_perf <- function(mem_perf) {
  mem_perf |>
    filter(mem_type %in% names(mem_types_report)) |>
    pivot_wider(
      names_from = mem_type,
      values_from = dprime
    ) |>
    column_to_rownames("subj_id") |>
    dist(method = "euclidean")
}

as_dist_vec <- function(vec, ..., size = NULL, diag = FALSE, upper = FALSE) {
  if (is.null(size)) {
    size <- 0.5 + sqrt(0.25 + 2 * length(vec))
  }
  stopifnot(all.equal(size, as.integer(size)))
  structure(
    vec,
    class = "dist",
    Size = size,
    Diag = diag,
    Upper = upper,
    ...
  )
}
