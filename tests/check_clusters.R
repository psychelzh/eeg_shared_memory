# Run from the repository root: Rscript tests/check_clusters.R
library(tidyverse)
source("R/clusters.R")

all_clusters <- find_cluster(c(2, 0, 2), c(0.01, 0.5, 0.01))
largest <- find_cluster(c(2, 0, 2), c(0.01, 0.5, 0.01), keep = "largest")
stopifnot(nrow(all_clusters) == 2L, nrow(largest) == 1L,
          largest$cluster_mass == 2)

unique_max <- find_cluster(c(2, 0, 3), c(0.01, 0.5, 0.01),
                           index = c(10, 20, 30), keep = "largest")
empty <- find_cluster(c(2, 0, 2), rep(0.5, 3), keep = "largest")
stopifnot(unique_max$cluster_mass == 3, unique_max$start == 30,
          nrow(empty) == 1L, empty$cluster_mass == 0, is.na(empty$start))

# Two permutations must contribute two maxima, even when one has tied clusters.
observed <- tibble(cca_id = 1L, time_id = 1:3,
                   statistic = c(3, 0, 0), p.value = c(0.01, 0.5, 0.5))
permuted <- tibble(cca_id = 1L, tar_rep = rep(1:2, each = 3),
                   time_id = rep(1:3, 2), statistic = c(2, 0, 2, 4, 0, 0),
                   p.value = c(0.01, 0.5, 0.01, 0.01, 0.5, 0.5))
result <- calc_clusters_stats(observed, permuted)
stopifnot(nrow(result) == 1L,
          identical(sort(result$cluster_mass_perm[[1]]), c(2, 4)),
          isTRUE(all.equal(result$p_perm, 2 / 3)))
cat("PASS: tied maxima, unique maxima, no clusters, and permutation p-value.\n")
