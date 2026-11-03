# EEG shared memory analysis

The current analysis is defined in `_targets.R`, with computation helpers in
`R/` and report plots in `R_quarto/`. `index.qmd` reads the computed targets.
Older exploratory scripts in `code/`, `_scripts/`, and `backup/` are research
history rather than alternative entry points for the current pipeline.

## Check and rerun

Use the project R environment and run commands from the repository root.
The behavioral checks require the local participant list and retrieval data.

```sh
Rscript tests/run_checks.R
```

This checks scoring, outcome selection, cluster permutation bookkeeping, and
the report's target references. It does not execute analysis targets or render
the report.

Run the analysis in R:

```r
targets::tar_make()
```

Workers are split into two crew pools. The `heavy` pool has one worker for
full-data EEG calculations and consumers of multi-GB dynamic pattern tables,
including their permutations, region construction, and reaction-time controls.
The `light` pool has 11 local workers (24 on SGE) for participant branches and
smaller statistics. The total budget remains 12 locally and 25 on SGE. Large
tasks can run alongside light work without blocking the main scheduler.
This limits heavy-task concurrency; it is not a machine-wide RAM reservation.

Cluster-permutation tasks now compute the statistic inside each permutation
replicate, before batching or writing its result. The cache retains statistics
rather than full permuted subject tables or Mantel model objects. ISPS
permutations similarly retain their means instead of subject-distance matrices.
Permutation counts, test formulas, and cluster correction remain unchanged.
Changed permutation target names change their deterministic RNG seeds: rerun
p-values can differ by Monte Carlo variation from historical outputs.

Large-data tasks use transient memory and pre-task garbage collection.
Permutation branches explicitly use `retrieval = "worker"`: workers read
their dependencies from the cache, avoiding the automatic main-process
retrieval/transport used for some dynamic branches. This is scoped to the
permutation factories and ISPS permutations, not a global policy change.
Local logs are under `logs/crew/light/` and `logs/crew/heavy/`. Default callr
isolation remains enabled. Full-pipeline speed and peak RAM still need to be
measured on the complete dataset; worker count alone does not guarantee speed.

An opt-in integration experiment uses tiny generated distances, two actual
worker pools and an isolated temporary cache. It checks fused Mantel results,
cluster results and overlapping execution, without running the main analysis:

```sh
Rscript tests/smoke_worker_pools.R
```

To isolate the repeatedly crashing GWS task without rerunning the full pipeline:

```r
targets::tar_make(names = data_gws_dynamic_alexnet_conv2d_5_11, use_crew = FALSE)
```

Worker crashes may not populate `.Last.error` in the interactive session.
Keep the printed task name and inspect the newest files under `logs/crew/`.
Do not increase crash retries as a substitute for diagnosing repeated failures.

Then regenerate the report and exported figures:

```sh
quarto render index.qmd
```

Quarto execution is not frozen: rendering reads the current target outputs.
Existing cached results and figures remain historical outputs until rerun.
There is no need to delete `_targets/` manually; `targets` tracks changed
commands and dependencies. Removed targets may leave unused cached objects.

## Scoring and statistics

- `calc_mem_perf()` computes Overall d-prime locally with the binary log-linear
  correction `(count + 0.5) / (N + 1)`.
- `calc_mem_perf_precise()` retains four-category smoothing with denominator
  `N + 2`, followed by the independent Remember/Know familiarity adjustment.
  Zero-frequency response categories contribute zero counts.
- `corr_mem()` defaults to the named `dprime` outcome. Multi-outcome tables
  require `col_perf`, such as `col_perf = dprime_rem`. Outcome rows must be
  unique per subject; unused forwarded arguments are errors.
- Each permutation contributes one maximum cluster mass, including tied
  maxima. Dynamic plot annotations use Holm-adjusted cluster p-values.

The retired mediation models and their disabled diagrams are no longer part
of the pipeline. Partial-representation correlations and the unique/shared
variance analysis remain in the report, including the existing aggregation
methods used for robustness comparisons. Historical mediation configuration
files remain in `config/`.
