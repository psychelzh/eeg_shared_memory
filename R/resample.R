resample <- function(n, size, paired = FALSE) {
  if (paired) {
    size <- size * 2
  }
  sampled <- sample.int(n, size)
  if (paired) {
    f <- rep(1:2, each = size / 2)
  } else {
    f <- rep(1, size)
  }
  split(sampled, f)
}
