# Constraint-v2 validation evidence

This directory contains compact release evidence, not a production scan.

- `assessed100.tsv`: deterministic 100-point zero-temperature reevaluation,
  with thermal results for all 22 eligible rows.
- `historical-verdicts.tsv`: old/new flags, numerical comparisons and per-row
  explanations. New vacuum/RG/thermal flags have no assumed historical values.
- `comparison-summary.json`: cross-host tolerances, counts, audit benchmarks,
  point-223615 crossing and integrated-pilot evidence.
- `amplitudes-{local,manto}.json`: 22 independently checked amplitude cases
  for each host, covering both micrOMEGAs versions.
- `runtime-{local,manto}.json`: 28 matched version/virtual-WZ benchmarks.
- `kernel-equivalence.json`: canonical source equivalence between installations.
- `tests-summary.json`: passing regression results and original log hashes.
- `stored100-{local,manto}.metadata.json`: unmodified provenance from the actual
  zero-temperature validation runs; source hashes retain their original meaning.
- `plots/v2-index.html`: all five new scientific plot families, PNG/PDF.

The input fixtures and deterministic selection provenance are in
`benchmarks/v2`. The 100-point sample is deliberately stratified; it is not a
statistical population estimate. All 100 points were reevaluated at zero
temperature on both hosts. The 22 eligible thermal points were run locally;
the explicit point-223615 and integrated-pilot thermal gates were also run on
manto. Validation Tmax is 300 GeV; the next-scan profile uses 1000 GeV.

Full native logs/traces and earlier validation attempts remain in the local
`validation/v2` tree and manto's isolated `runtime-v2/source/validation/v2` tree.
The committed artifacts omit large generated campaign products. Runtime
receipts live outside Git and record the final deployed source commit.

Reproduce the comparisons after collecting the native outputs using
`tools/summarize_validation_v2.py`. See `docs/next-scan-readiness-v2.md` for the
validation interpretation and launch commands.
