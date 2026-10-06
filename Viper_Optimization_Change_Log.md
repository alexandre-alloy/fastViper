# fastViper optimization change log

This log summarizes local changes to the viper 1.46.0 source package.

## Performance and enrichment

- Replaced row and column ranking `apply()` calls with `matrixStats` rank functions, retaining average tie handling.
- Replaced the full-matrix imputation and multiplication used by `frvarna()` with compiled row variance calculation while preserving missing-value behavior.
- Reduced allocations in the common complete-data aREA loop by avoiding an all-ones weight matrix and redundant mask multiplications.
- Added `enrichment.method = c("aREA", "NaRnEA")` to `viper()`. aREA remains the default. NaRnEA returns NES normally and PES when `nes = FALSE`; it rejects missing signature values, as in the upstream algorithm. With NaRnEA, `dnull` is not used for primary score normalization.
- Added a multi-sample NaRnEA implementation using sparse regulator matrices and shared ranking/covariance calculations.
- Vectorized the common complete-data, non-robust lineal/rank `viperRPT()` calculation; unsupported or degenerate inputs retain the original fallback.
- Cached directed target intersections in `shadowRegulon()`.

## Multicore execution

- Added a common parallel helper: PSOCK on Windows and fork workers on Unix-like systems. Existing `cores` arguments activate it when greater than one; serial remains the default.
- Applied it to VIPER, aREA, NaRnEA, signatures, bootstrap/null-model routines, shadow/pleiotropy routines, msVIPER workflows, and `viperRPT()`.
- Added a nested-call guard so only one worker layer is active at a time. For example, when VIPER distributes regulons across workers, aREA runs serially within each worker.
- On Windows, workers receive a reduced closure and the explicit inputs needed for their tasks. They still hold serialized process-local copies; matrices are not shared or memory mapped.
- Added an estimated aggregate input-copy budget controlled by `options(viper.max.parallel.memory = ...)`, defaulting to 64 GiB. This estimate does not cap total process memory.
- Added minimal dependency capture so the helper also works when code is loaded with `source()`.

## Validation performed

- Compared optimized calculations with original viper 1.46.0 on targeted workloads. The 50-patient signature, null model, and activity outputs matched within 1e-10; the signature calculation improved from 7.67 s to 2.90 s.
- A 1,080-patient serial workflow completed with 100 permutations in `viperSignature()`; this full workflow has not been benchmarked with multiple workers.
- NaRnEA scores matched upstream matrixNaRnEA within floating-point tolerance on synthetic tests. A 10,000-gene, 40-sample test was 2.39x faster.
- Windows PSOCK checks matched serial aREA, VIPER, metaVIPER, NaRnEA, and signature outputs within 1e-12. Repeated seeded signature calculations were identical. Installed-package and `source()` helper smoke checks passed.
- The common aREA optimization was 1.52x faster in a 30,000-gene, 100-sample synthetic benchmark. The vectorized lineal/rank RPT paths were substantially faster on a 600-gene, 80-sample synthetic benchmark and agreed within 1e-10.
