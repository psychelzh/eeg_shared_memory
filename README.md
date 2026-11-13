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
the report's target references. A tiny temporary pipeline checks cache relocation;
the checks do not execute the study analysis or render the report.

After a completed analysis, refresh the report and manuscript statistics without
rebuilding the targets pipeline:

```sh
quarto render index.qmd --to html --execute
Rscript scripts/summarize_completed_results.R
```

The report refreshes `figures/`, `_freeze/index/` and `_site/index.html`.
The summary script writes numeric tables to a local, Git-ignored
`docs/analysis-update-YYYY-MM-DD/` directory. Pass a directory as its first
argument to choose another output location. Generated exports and personal
manuscript notes are not part of the repository.

Run the analysis in R:

```r
targets::tar_make()
```

The pipeline uses shared native partition targets and a common task-building framework.
Pattern tables are materialized by CCA component or brain region and reused by
IGS, ISS/IWS, ISPS, GSS/GWS, SMC/partial Mantel, region and RT analyses. A single
worker prepares keys and native pieces; calculations share 12 local workers (32 on
SGE). Total worker budgets are 13 locally and 33 on SGE. The local default is
tuned for this workstation's i7-12700F (12 cores / 20 threads) and 64 GiB RAM.
Local workers inherit single-thread limits for OpenMP/Arrow, OpenBLAS, MKL and
RcppParallel to avoid nested CPU parallelism. Mantel already uses one process.
The earlier Mantel benchmark favoured 14, but a subsequent full-run startup
failure prompted a more conservative default of 12. Local workers have a
120-second launch allowance and inherit `RENV_CONFIG_SANDBOX_ENABLED=FALSE`:
one failed run accumulated hundreds of workers while renv sandbox activation
took about 2,900 seconds. Project renv activation and its library remain in use;
the extra system-library sandbox is skipped for new workers. This is not a
full-pipeline memory or speed guarantee. To change the local compute pool,
use `Sys.setenv(EEG_WORKERS = "12")`.
Changing these settings does not resize an already-running controller.
Runtime setup is kept in `R/tar_runtime.R`; `_targets.R` calls
`tar_source()` once to load the R functions, then `configure_targets_runtime()`.

SGE uses serial one-slot jobs, with 6 GB requested per compute worker and 24 GB
for preparation. It reads `qconf -sc` before submission and selects a requestable
`m_mem_free` or `mem_free` resource. Set `EEG_SGE_MEMORY_RESOURCE` for a custom
resource. The supplied cluster configuration has nonconsumable `mem_free`, so
these are free-memory thresholds, not reservations. In that case jobs also
require `mem_total=400G`, selecting the `cogmap`/`schema` class of hosts from the
supplied qhost snapshot. Override that floor with `EEG_SGE_MIN_HOST_GB` (0 disables
it). Missing resources stop configuration before submission. `EEG_SGE_WORKERS`
overrides the compute count
(default 32). Job scripts explicitly export the same single-thread and renv
sandbox settings, with per-job logs under `logs/crew/sge/light/` and
`logs/crew/sge/heavy/`. Idle workers are retained for up to 300 seconds to reuse
their loaded environment. Actual concurrency depends on queue access and quotas.

Computation workers retrieve only their persisted native partitions, once per batch,
and compute complete-subject
time windows sequentially within that job. Parallelism is across partitions
and batches. The cache retains permutation statistics or ISPS means. Full time
series are combined before cluster inference. Raw EEG regions and full CCA
pattern construction are now separate computation branches. RT memory-label
permutations reuse cached numeric IGS/ISS tables. Small numeric tables keep
whole-table batch parallelism; joint prediction models keep all components.

The shared seed plan preserves cross-window memory/target shuffles and
independent subject/window pattern shuffles. Draws are stable across batching
and worker assignment. This establishes new Monte Carlo draws relative to
historical caches, so affected p-values can change. Outer counts remain 1,000;
Mantel inner counts remain 9,999 with no nested process parallelism. Mantel
batches use two outer replicates, reducing the work lost on a worker crash.

Read the [architecture and reuse guide](docs/partitioned-pipeline.md) for the
registry, builders, randomisation rules and resource policy. Only preparation
uses the single-worker pool; this is not a machine-wide RAM reservation.
Native partition payloads live under `_targets/objects/`, with no manifest,
generation directory or machine path in jobs. The keys and each piece are prepared
serially in the heavy pool. Preparing a piece may reread an entire upstream table
or EEG input; this is bounded preparation work, not repeated per permutation batch.
Relocating a complete native cache preserves identities with compatible inputs
and environments. Upgrading the completed file-backed cache changes dependencies
and can require recomputation; it is not an automatic metadata migration.
See the guide's project relocation section before running the completed cache.
No full-pipeline speed or peak-memory guarantee is implied by the small tests.

The following opt-in checks use isolated temporary stores:

```sh
Rscript tests/smoke_partitioned_pipeline.R
Rscript tests/smoke_real_partitions.R
Rscript tests/smoke_native_partitions.R
```

The first uses generated inputs and checks actual workers, file-to-native upgrade,
reused native objects, mapped outcomes and all major calculation families. The second
requires cached project data and checks actual GWS distances with 9,999 inner
permutations. The earlier [ISPS pilot](docs/isps-chunk-experiment.md) records the
initial design and its benchmark limits. The native-partition prototype checks
per-CCA object storage, worker execution, numerical equivalence, relocation and
native pruning. The [pipeline guide](docs/partitioned-pipeline.md) describes the
production mechanism and its validation limits.

Default callr isolation remains enabled. Local worker logs are under
`logs/crew/light/` and `logs/crew/heavy/`. To run one analysis first:

```r
targets::tar_make(names = clusters_stats_isps_dynamic)
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
