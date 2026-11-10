# Run from the repository root: Rscript tests/run_checks.R
checks <- list.files("tests", pattern = "^check_.*\\.R$", full.names = TRUE)
for (check_file in checks) {
  cat("Checking", check_file, "\n")
  sys.source(check_file, envir = new.env(parent = globalenv()))
}
cat(
  "All",
  length(checks),
  "checks passed. No study analysis targets were executed.\n"
)
