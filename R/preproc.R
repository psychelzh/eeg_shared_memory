# keep pairs of the same response to check the inter-subject similarity
extract_response_shared <- function(events_encoding, events_retrieval) {
  events_retrieval |>
    filter(memory_type != 0) |>
    pivot_wider(
      id_cols = word_id,
      names_from = subj_id,
      values_from = response_type
    ) |>
    mutate(
      dplyover::across2x(
        -1, -1,
        ~ if_else(.x == .y, .x, NA),
        .comb = "minimal"
      )
    ) |>
    select(contains("_")) |>
    pivot_longer(
      -word_id,
      names_to = "subj_pair",
      values_to = "response_type_shared"
    ) |>
    separate(
      subj_pair,
      c("subj_id_col", "subj_id_row"),
      convert = TRUE
    ) |>
    mutate(
      response_type_shared = case_match(
        response_type_shared,
        "remember" ~ "Rem",
        "know" ~ "Know",
        "unsure" ~ "Unsure",
        "new" ~ "New",
        .ptype = factor(levels = c("Rem", "Know", "Unsure", "New"))
      )
    ) |>
    nest(.by = word_id, .key = "resp_matched") |>
    left_join(
      distinct(events_encoding, trial_id, word_id, word_category),
      by = "word_id"
    )
}

filter_shared <- function(file, response_shared) {
  arrow::read_parquet(file) |>
    inner_join(response_shared, by = "trial_id") |>
    mutate(
      filtered = map2(
        resp_matched,
        fisher_z,
        ~ .x |>
          add_column(fisher_z = .y) |>
          filter(!is.na(response_type_shared))
      ),
      .keep = "unused"
    ) |>
    unnest(filtered) |>
    filter(
      sum(!is.na(fisher_z)) >= 5,
      .by = c(region_id, word_category, contains("subj_id"))
    )
}

# average inter-subject similarity across trials for trial-level analysis
average_rs_trials <- function(file_parquet,
                              col_rs = fisher_z,
                              col_trial = trial_id,
                              scalar_rs = FALSE) {
  dat <- arrow::read_parquet(file_parquet) |>
    nest(.by = -c({{ col_trial }}, {{ col_rs }}))
  if (scalar_rs) {
    dat |>
      mutate(
        mean_fisher_z = map_dbl(
          data,
          ~ mean(pull(.x, {{ col_rs }}), na.rm = TRUE)
        ),
        .keep = "unused"
      )
  } else {
    dat |>
      mutate(
        mean_fisher_z = map(
          data,
          ~ do.call(rbind, pull(.x, {{ col_rs }})) |>
            colMeans(na.rm = TRUE)
        ),
        .keep = "unused"
      )
  }
}


# prepare shuffled behavioral measures by permuting subject id
permutate_behav <- function(data, cols_id) {
  data_ids <- unique(data[cols_id])
  data_ids_perm <- data_ids[sample.int(nrow(data_ids)), ]
  suff_tmp <- "_perm"
  names(data_ids_perm) <- paste0(cols_id, suff_tmp)
  bind_cols(data_ids, data_ids_perm) |>
    left_join(data, by = cols_id) |>
    select(-all_of(cols_id)) |>
    rename_with(
      ~ str_remove(.x, suff_tmp),
      ends_with(suff_tmp)
    )
}

permutate_simil <- function(simil) {
  perm <- sample.int(attr(simil, "Size"))
  as.dist(as.matrix(simil)[perm, perm])
}
