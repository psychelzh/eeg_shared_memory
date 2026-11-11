# Graph-building configuration, loaded by tar_source(). Runtime target commands
# depend on the selected native branches, not on this registry as a whole.
partition_sources <- local({
  sources <- c(
    "patterns_indiv_dynamic",
    "patterns_group_dynamic_loo",
    "patterns_indiv_whole",
    "patterns_group_whole_loo",
    "patterns_group_dynamic",
    "patterns_group_whole",
    "patterns_indiv_dynamic_regions",
    "patterns_indiv_whole_regions",
    "patterns_indiv_dynamic_rt_locked",
    "patterns_indiv_whole_rt_locked",
    "patterns_group_whole_loo_rt_locked",
    "data_isps_dynamic",
    "data_isps_whole",
    "data_isps_partial_group_dynamic",
    "data_isps_partial_group_whole",
    "data_isps_partial_semantic_dynamic",
    "data_isps_partial_semantic_whole"
  )
  stats::setNames(paste0("partitions_", sources), sources)
})
