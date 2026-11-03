# Run from the repository root: Rscript tests/check_corr_mem.R
library(dplyr)
source("R/prediction.R")

# Opposite correlations expose accidental selection of the last outcome column.
predictors <- tibble(subj_id = 1:6, iss = 1:6)
outcomes <- tibble(
  subj_id = 1:6,
  dprime_rem = c(6, 4, 5, 2, 3, 1),
  dprime_know = c(1, 3, 2, 5, 4, 6)
)
for (outcome in c("dprime_rem", "dprime_know")) {
  actual <- corr_mem(
    predictors, outcomes, col_perf = all_of(outcome), alternative = "greater"
  )
  expected <- cor.test(
    predictors$iss, outcomes[[outcome]], alternative = "greater"
  )
  stopifnot(
    isTRUE(all.equal(unname(actual$estimate), unname(expected$estimate))),
    isTRUE(all.equal(actual$p.value, expected$p.value)),
    actual$alternative == "greater"
  )
}

# Unknown arguments must fail even when forwarded to cor.test() via dots.
for (argument in c("name_perf", "alernative")) {
  args <- list(data = predictors, mem_perf = outcomes)
  args[[argument]] <- "dprime_rem"
  error <- tryCatch(do.call(corr_mem, args), error = identity)
  stopifnot(inherits(error, "error"), grepl(argument, conditionMessage(error)))
}

# Overall memory is explicit by default, regardless of column order.
overall <- mutate(outcomes, dprime = dprime_rem, .after = subj_id)
actual <- corr_mem(predictors, overall, method = "spearman", exact = FALSE)
expected <- cor.test(predictors$iss, overall$dprime,
                     method = "spearman", exact = FALSE)
stopifnot(isTRUE(all.equal(unname(actual$estimate), unname(expected$estimate))))

# Multi-outcome data require an explicit outcome; joins must not duplicate subjects.
missing_outcome <- tryCatch(corr_mem(predictors, outcomes), error = identity)
duplicate_outcomes <- bind_rows(overall, slice(overall, 1))
duplicate_error <- tryCatch(corr_mem(predictors, duplicate_outcomes), error = identity)
multiple_error <- tryCatch(
  corr_mem(predictors, outcomes, col_perf = starts_with("dprime")), error = identity
)
stopifnot(inherits(missing_outcome, "error"), inherits(duplicate_error, "error"),
          inherits(multiple_error, "error"))

# Explicit selection continues to forward supported cor.test arguments.
actual <- corr_mem(predictors, outcomes, col_perf = dprime_know,
                   method = "spearman", exact = FALSE)
expected <- cor.test(predictors$iss, outcomes$dprime_know,
                     method = "spearman", exact = FALSE)
stopifnot(isTRUE(all.equal(unname(actual$estimate), unname(expected$estimate))))
cat("PASS: outcome selection, forwarded arguments, and unused-dot errors.\n")
