source("R/tar_runtime.R")

complexes <- c(
  "#name shortcut type relop requestable consumable default urgency",
  "mem_free mf MEMORY <= YES YES 0 0",
  "m_mem_free mmf MEMORY <= YES YES 0 0",
  "mem_total mt MEMORY <= YES NO 0 0",
  "site_memory sm MEMORY <= FORCED JOB 0 0",
  "slots s INT <= YES YES 1 1000"
)
stopifnot(
  select_sge_memory_config(complexes)$resource == "m_mem_free",
  select_sge_memory_config(complexes, "mem_free")$resource == "mem_free",
  select_sge_memory_config(complexes, "site_memory")$resource == "site_memory",
  select_sge_memory_config(complexes[-3L])$resource == "mem_free",
  select_sge_memory_config(complexes)$consumable,
  select_sge_memory_config(complexes)$min_host_gigabytes == 0,
  sge_memory_directive("mem_free", 6) == "#$ -l mem_free=6G",
  sge_memory_directive("m_mem_free", 24) == "#$ -l m_mem_free=24G"
)
for (resource in c("slots", "missing")) {
  error <- tryCatch(
    select_sge_memory_config(complexes, resource),
    error = identity
  )
  stopifnot(inherits(error, "error"))
}
nonconsumable <- "mem_free mf MEMORY <= YES NO 0 0"
site <- c(
  nonconsumable,
  "mem_total mt MEMORY <= YES NO 0 0",
  "slots s INT <= YES YES 1 1000"
)
config <- select_sge_memory_config(site)
stopifnot(
  config$resource == "mem_free",
  !config$consumable,
  config$min_host_gigabytes == 400,
  select_sge_memory_config(site, min_host_gigabytes = 0)$min_host_gigabytes ==
    0,
  sge_memory_directive("mem_total", config$min_host_gigabytes) ==
    "#$ -l mem_total=400G"
)
for (invalid in list(-1, NA_real_, Inf)) {
  error <- tryCatch(
    select_sge_memory_config(site, min_host_gigabytes = invalid),
    error = identity
  )
  stopifnot(inherits(error, "error"))
}
for (lines in list(character(), "# header only", nonconsumable)) {
  error <- tryCatch(select_sge_memory_config(lines), error = identity)
  stopifnot(inherits(error, "error"))
}
cat(
  "PASS: SGE resource selection, real-site nonconsumable fallback, host-memory filter and 6G/24G directives.\n"
)
