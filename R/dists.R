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

