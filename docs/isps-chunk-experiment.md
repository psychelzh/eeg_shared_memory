# Dynamic ISPS component experiment

This records the initial ISPS-only pilot. The current implementation generalises
the approach across analysis families; see [the shared framework](partitioned-pipeline.md).
The numbers below describe that initial pilot rather than a new full-run estimate.

Local experiment on 2026-11-03, using the project's restored R environment.
The main analysis was stopped. All experiments used temporary files/stores;
the formal analysis cache was read as input and not recomputed by the experiments.

The old calculation constructs all permuted subject patterns before calculating
ISPS. The new calculation reads one CCA component and generates the trial-label
orders in original row order, then computes one time window at a time. Component
and replicate seeds are stored explicitly, with fixed RNG settings. Subjects
remain together within each window, and complete time series feed cluster tests.

## Real-data pilot

206 subjects, three CCA components, 12 time windows, one permutation per component.
Both implementations used identical permutation orders. The reference ran in
one process; the new implementation used three concurrent crew workers.
Component files were prepared before timing either calculation.

| Measurement | Reference | Component workers |
|---|---:|---:|
| Calculation interval | 33.74 s | 12.86 s |
| Process/pipeline time, including startup and input loading | 37.00 s | 37.61 s |
| Largest individual process peak working set | 2.21 GiB | 0.62 GiB |

The calculation interval was about 2.6 times faster. A single-permutation pilot
is too small to amortise worker startup, so its total time was similar.
The per-process memory figures are not total machine RAM: all three worker
peaks sum to approximately 1.85 GiB, excluding the scheduler/driver.
All 36 permutation-summary rows matched the reference.

A separate complete-component check retained all 206 subjects and 52 windows:
one permutation took 49.78 seconds including file loading, with a 1.65 GiB
process peak working set. Observed distances exactly matched a fresh evaluation
of the original function. Historical cached distances differed by at most
5.6e-17 numerically and retained the old `use = "pairwise"` call attribute.

## Validation and scope

The current `tests/check_partitioning.R` checks observed distances, missing values, fixed
permutation equivalence, invariance to batch boundaries, null ranks, pointwise
p-values and time-cluster results. `tests/smoke_partitioned_pipeline.R` upgrades a cached
legacy pipeline in an isolated store and checks file tracking, six cross
branches, overlapping workers and all permutation results. It also verifies
that deleting a component file invalidates its producer.

The ISPS pilot uses three component files and 100 batches, creating 300 jobs with
1,000 permutations per component/window. These jobs share the existing light
pool; eleven simultaneous ISPS workers and the full 1,000-permutation run were
not benchmarked. Mantel's nested permutation counts are retained; the shared
framework now also partitions other large-table analyses. This experiment does not
establish a new completion estimate for the entire pipeline.
