# Native targets partitions and cache cleanup

## Assessment

The production pipeline now replaces custom partition files with persistent
native dynamic targets, using the factories in `R/tar_partitions.R`.
Each source partition should be a real branch whose value is the complete table
for one CCA component or brain region. Downstream computation can cross these
branches with the recorded permutation batches. `targets` then owns the stored
objects, metadata, dependency tracking and ordinary pruning, without a custom
manifest or `generation-*` directory in the generated target values.

The official [dynamic branching guide](https://books.ropensci.org/targets/dynamic.html)
describes dynamic branches and grouped iteration. The
[target options reference](https://docs.ropensci.org/targets/reference/tar_target.html)
explains worker retrieval and static upstream objects versus persisted branches.

## Executed prototype

```sh
Rscript tests/smoke_native_partitions.R
```

The opt-in prototype uses generated data in temporary projects and two local
crew workers. It does not run or migrate the completed study pipeline.

- Three CCA source partitions persist as three independent objects under
  `_targets/objects/`, with list iteration and worker retrieval.
- Three component branches crossed with four batches produce 12 computation
  branches. Each receives 32 source rows for its own component, rather than the
  combined 96-row table. Each window still contains all eight synthetic subjects.
- Eight recorded outer replicates produce 96 statistic rows for each of ISS
  and IGS. Both families match the current file-backed implementation to
  tolerance 1e-12, using the same seeds and randomisation helpers.
- Copying the completed project/cache to another root reports zero outdated
  targets and reuses identical saved results.
- Removing an obsolete target and calling `tar_prune(cloud = FALSE)` removes its
  native object and metadata while preserving the retained results.
- No partition files are generated under `_targets/user/` by the native design.

These checks passed. They demonstrate the mechanism, not a full-study RAM or
throughput benchmark, a real Windows-to-Linux transfer, or equivalence testing
for every Mantel/ISPS/region analysis.

## Production implementation and validation

Many existing calculations already branch by CCA and can reuse those native
branches. Sources that branch by subject or return whole tables need preparation
that materializes complete CCA/region subsets before the permutation batches.

Simply marking a large static whole table as `iteration = "group"` does not
establish efficient disk partitioning. With worker retrieval, a static upstream
object can be reread in full for many branches. Materialize grouped subsets once,
then branch the expensive computations over those persisted pieces. Preparation
must retain the one-worker heavy resource policy. Manager memory and repeated
whole-table reads need measurement for sources not already divided appropriately.

Multiple sources must be aligned by validated CCA/region keys. `map()` aligns
branch positions rather than performing a relational join on these keys. The
prototype checks component identity explicitly. Preserve complete-subject windows,
shared outer seeds and keyed row permutations.

The production registry and raw EEG region preparation now use this mechanism.
There are 842 target definitions, including 18 key catalogs and 18 native
partition roots. None of the production partition targets uses `format = "file"`.
Jobs contain keys only. Source catalogs must have identical ordered keys, and
each piece is checked against its job before computation. Existing seed plans,
batch sizes, complete-subject windows and heavy/compute resource policies remain.

Validation uses temporary stores and leaves the completed study cache untouched:

- `tests/check_native_partitions.R` tests the actual production factories,
  relocation without recomputation, missing-object repair, real-change
  invalidation and key-alignment failures.
- `tests/smoke_partitioned_pipeline.R` upgrades a file-backed fixture to native
  partitions and compares its saved ISS null results. It also checks parallel
  ISS/IGS/ISPS/Mantel, partial Mantel, mapped outcomes, six raw EEG regions,
  overlapping workers and missing-object invalidation.
- `tests/smoke_real_partitions.R` checks native pieces against monolithic
  calculations on cached real study data, including 9,999 inner Mantel
  permutations.

These checks passed. They are not a full-study memory/throughput benchmark or an
actual Windows-to-Linux transfer. Preparation still reads a whole source table
for each key, or the raw EEG array for each region, with one heavy worker.
That increases preparation I/O; subsequent permutation batches retrieve only
their persisted piece.

This is a pipeline/cache migration: dependencies and branch identities change.
No full study run or automatic old-cache conversion has been performed. The next
full run can rebuild affected expensive targets. Legacy file helpers remain for
old-cache/reference tests; existing referenced files must be retained until the
native run has completed and its results have been checked. The cleanup counts
below describe the completed file-backed graph before this upgrade.

## Completed local cleanup

The audit retained the current graph's recorded children, referenced generated
file outputs, and file references embedded in current jobs. Copied HPC paths
were mapped to local counterparts for verification; saved results were unchanged.

| Item | Before | After |
|---|---:|---:|
| Cached object files | 51,685 | 33,034 |
| Cached object size | 231.737 GiB | 18.318 GiB |
| Generated partition files | 150 | 81 |
| Generated partition size | 28.625 GiB | 14.313 GiB |
| D drive free space | 125.563 GiB | 353.339 GiB |

The audit removed 18,720 files with 227.731 GiB total logical size, then pruned
obsolete metadata with the native API. All retained files passed existence and
size checks. All 221 checked final-result hashes are unchanged. All 824 current
definitions retain metadata. All 17 final GWS/GSS/SMC cluster outputs can be read
and retain 1,000 null maxima per row.

Some empty generation directories remain: automatic approval review blocked the
additional empty-directory removal command. This does not retain deleted payload
or affect verified space recovery.

The deletion inventory, fingerprints, metadata backup and completion record are
kept locally in `logs/cache-cleanup-20261110/` (ignored by Git). The cleanup script
is specific to that audited inventory and checks source/metadata fingerprints.
