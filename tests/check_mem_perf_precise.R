# Run from the repository root: Rscript tests/check_mem_perf_precise.R
library(dplyr)
library(tidyr)
source("R/behav.R")

# Independent reference: build the full contingency table explicitly.
reference <- function(data) {
  bind_rows(lapply(sort(unique(data$subj_id)), function(id) {
    x <- data[data$subj_id == id, ]
    counts <- table(
      factor(x$old_new, levels = 1:2),
      factor(x$resp, levels = 1:4)
    )
    rates <- (counts + 0.5) / (rowSums(counts) + 2)
    rem <- rates[, 1]
    know <- rates[, 2] / (1 - rem)
    tibble(
      subj_id = id,
      dprime_rem = qnorm(rem[1]) - qnorm(rem[2]),
      dprime_know = qnorm(know[1]) - qnorm(know[2])
    )
  }))
}

check <- function(data) {
  actual <- arrange(calc_mem_perf_precise(data), subj_id)
  expected <- reference(data)
  stopifnot(
    nrow(actual) == length(unique(data$subj_id)),
    all(is.finite(as.matrix(actual[-1]))),
    isTRUE(all.equal(as.data.frame(actual), as.data.frame(expected),
                     check.attributes = FALSE, tolerance = 1e-12))
  )
  actual
}

# Missing Remember false alarms; no Know responses at all; extreme responses.
synthetic <- tibble(
  subj_id = rep(1:3, each = 8),
  old_new = rep(rep(1:2, each = 4), 3),
  resp = c(1, 2, 3, 4, 2, 2, 3, 4,
           1, 1, 3, 4, 1, 3, 4, 4,
           1, 1, 1, 1, 4, 4, 4, 4)
)
invisible(check(synthetic))
invisible(check(synthetic[nrow(synthetic):1, ]))

events <- read.delim("data/behav/retrieval.tsv")
ids <- scan("data/subj_206.txt", quiet = TRUE)
events <- events[events$subj %in% ids & events$resp != 0, ]
events$subj_id <- match(events$subj, ids)
actual <- check(events)
stopifnot(nrow(actual) == 206L)

# Where all four response cells existed, the correction must be unchanged.
complete_ids <- events |>
  count(subj_id, old_new, resp) |>
  count(subj_id, name = "cells") |>
  filter(cells == 8L) |>
  pull(subj_id)
invisible(check(filter(events, subj_id %in% complete_ids)))

events$score <- as.double(xor(events$old_new == 1, events$resp >= 3))
combined <- left_join(calc_mem_perf(events), actual, by = "subj_id")
print(cor(as.data.frame(combined[-1])))
cat("PASS: zero-count cases, reversed rows, and all 206 participants.\n")
