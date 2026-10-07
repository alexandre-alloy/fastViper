# fastViper

`fastViper` is a locally optimized version of the R `viper` package for inferring regulator activity from gene expression data. It keeps the package's existing VIPER and msVIPER workflows and adds performance improvements, an optional NaRnEA enrichment method, and cross-platform multicore execution. The package source directory remains named `Viper_Optimized`; `fastViper` is the project/repository name.

The package is based on viper 1.46.0. See [`Viper_Optimization_Change_Log.md`](Viper_Optimization_Change_Log.md) for the changes and validation results.

## What it does

- **VIPER and msVIPER:** infers transcription factor or regulator activity from an expression matrix and a regulon/network.
- **aREA enrichment:** the default enrichment method, with faster ranking and lower-overhead paths for common complete-data workloads.
- **NaRnEA enrichment:** an optional analytic enrichment method that can be selected in `viper()` with `enrichment.method = "NaRnEA"`. With `nes = FALSE`, it returns NaRnEA's proportional enrichment scores (PES); otherwise it returns normalized enrichment scores (NES).
- **Expression signatures:** `viperSignature()` creates signatures against a reference, including z-score, t-test, and mean modes.
- **Other optimized paths:** rank and variance calculations use `matrixStats`; a common non-robust `viperRPT()` path uses vectorized weighted least squares; `shadowRegulon()` scores directed regulator pairs from cached overlap indices and can distribute focal regulators across workers.
- **Optional multicore:** `cores = 1` runs serially. Values above one enable PSOCK workers on Windows and fork workers on Unix-like systems for supported independent tasks.

## Install

Install from a local checkout with R:

```r
remotes::install_github("alexandre-alloy/fastViper")
library(fastViper)
```

## Basic use

```r
library(fastViper)

# expr: genes in rows, samples in columns
# regulon: a regulon object or a list of regulons
activity <- viper(
  eset = expr,
  regulon = regulon,
  method = "none",
  enrichment.method = "aREA"
)

# Use analytic NaRnEA enrichment instead
activity_narnea <- viper(
  eset = expr,
  regulon = regulon,
  method = "none",
  enrichment.method = "NaRnEA"
)

# Enable multiple workers (including on Windows)
activity_parallel <- viper(
  eset = expr,
  regulon = regulon,
  method = "none",
  cores = 4
)
```

To compute sample signatures against a reference matrix:

```r
signatures <- viperSignature(
  eset = test_expr,
  ref = reference_expr,
  method = "zscore",
  per = 100,
  cores = 4
)
```

The resulting signatures can then be passed to `viper()` as the expression input. Choose `method` for `viper()` independently; for example, use `method = "none"` when the supplied input already contains the desired signatures.

## Multicore behavior and memory

The `cores` argument is available in `viper()` and several computationally intensive functions, including `aREA()`, `viperSignature()`, bootstrap/null-model routines, msVIPER workflows, and `viperRPT()`. Work is split across independent samples, regulons, or permutations where supported.

Nested calls use one worker layer: when an outer operation has parallelized its work, inner operations such as `viper()` calling `aREA()` run serially within each worker. This avoids multiplying worker pools and their memory use.

Windows PSOCK workers are separate R processes. They need serialized copies of their input data; this implementation does not use shared-memory or memory-mapped matrices. A worker limit estimates aggregate input-copy memory and defaults to 64 GiB. Adjust it for your machine with:

```r
options(viper.max.parallel.memory = 384 * 1024^3) # 384 GiB example
```

This is an estimate of input copies, not a cap on total R process memory. Outputs, temporary calculations, R overhead, and other applications also use RAM. If the requested worker count exceeds the estimated budget, the number of workers is reduced; if necessary the call runs serially.

## Validation

The local change log records comparisons against the unmodified viper 1.46.0 source, including matching outputs for optimized calculations, NaRnEA comparisons against the upstream implementation, and Windows PSOCK checks against serial results. A 30-patient, 100-regulator `shadowRegulon()` benchmark matched exactly, with a 1.18x serial speedup and 1.96x speedup using four workers. A larger 120-patient, 187-regulator run was 12.75x faster with 32 workers, again with identical outputs; sampled aggregate R-process memory peaked at 13.71 GiB. These checks cover targeted synthetic workloads and selected package workflows. A complete 1,080-patient analysis was run serially for performance evaluation; that full workflow has not been benchmarked with multiple workers.

## License and attribution

This project is based on the viper package. The original license and package notices are included in [`LICENSE`](LICENSE). Please retain those notices and review the license terms before redistribution.
