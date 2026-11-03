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
