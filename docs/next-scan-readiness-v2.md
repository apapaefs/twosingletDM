# Next-scan readiness: constraint pipeline v2

## Release scope

The active `vx=0` pipeline implements the approved plan. The next-scan profile
uses micrOMEGAs **7.1.4**, the verified observed **LZ WS2024** SI numerical table,
and a common **BSMPT 3.2.1** source on the laptop and manto. The corrected
micrOMEGAs 6.1.15 backend remains available for comparisons. Every evaluated
point is retained, including DM exclusions and unavailable assessments.

The release commit is the repository `HEAD` recorded in each site's external
`runtime-v2/runtime-manifest.json`. That receipt is generated after committing
and deploying; it contains executable hashes, compiled source hashes, package
versions, effective inputs and the release source commit. Committed validation
metadata keeps its original pre-release revision and individual source hashes.
It is never rewritten to imply that a run used a later commit.
Manto's staging export had no Git checkout; its recorded file hashes identify
the validated source before deployment into the production checkout.

The [pipeline specification](constraints-v2.md) defines the approximations,
field thresholds, status semantics, fingerprints and legacy handling. The
[configuration](../config/next-scan-v2.json) specifies 10,000 draws over four
seeds and a 2 GeV annihilation-resonance window. The full scan is not part of
release validation.

## Validation gates

| Gate | Evidence |
|---|---|
| Regression suite | 343 tests pass locally and on manto; no skipped archived-data regression |
| Theory counterexamples | Deeper X-breaking vacuum, finite but nonunitary running, repeated eigenvalues, degenerate/flat branches |
| Independent DM amplitudes | 22 cases per host, covering both backends, poles, mixing, interference, thresholds, on-shell widths and raw CalcHEP AA/GG normalization |
| Backend/runtime comparison | 28 matched version/virtual-WZ cases on both hosts; nonzero correction-hook counts in relic and indirect calculations |
| Stored points | 100 deterministic points reevaluated at zero temperature on each host; all comparable verdicts agree |
| Thermal validation | All 22 eligible stored points processed locally at Tmax = 300 GeV; point 223615 also evaluated independently on manto |
| Integrated pilot | Interruption/resume reproduces uninterrupted generation; changed physics settings reject resume before appending; a DM-failing point still reaches BSMPT and retains its GW flag |
| Concurrent campaign | Two isolated workers, four retained input rows, actual output manifests and candidate summaries |
| Plot production | All five new families, PNG/PDF and HTML index; the comprehensive suite also produced 69 figures |

The detailed [comparison summary](../validation/v2/release/comparison-summary.json)
contains audit results, tolerances and pilot evidence. The independent
[amplitude report](../validation/v2/release/amplitudes-local.json),
[runtime benchmarks](../validation/v2/release/runtime-local.json),
[100-point ledger](../validation/v2/release/assessed100.tsv), and
[old/new verdict table](../validation/v2/release/historical-verdicts.tsv) are
committed alongside the input fixtures.

### Numerical agreement and tolerances

For all 100 rows, local/manto values of relic abundance, Xf, Tf, SI cross section,
backend mediator widths and HiggsSignals chi-square agree exactly at the saved
precision. The 28 matched runtime cases also agree at the saved precision.
The acceptance tolerances are relative `2e-6` for the DM quantities (absolute
`1e-30`), and relative `1e-8`, absolute `1e-10` for HiggsSignals. Thermal
comparisons allow 0.003 GeV for transition temperatures and relative `1e-3`
for jump ratios. Refined crossings must have a bracket no wider than 0.001 GeV.

Independent loop checks use relative `2e-8`; on-shell partial-width comparisons
allow `2e-6` because the upstream branching-ratio text is rounded. The maximum
7.1.4-versus-6.1.15 relic difference across all tested settings is **0.125%**
(h1 pole, historical virtual-WZ-off setting); with the new virtual-WZ-on profile,
the largest difference is **0.0817%**. These are benchmark observations, not a
global uncertainty estimate. At the W-threshold benchmark, enabling virtual
W/Z channels changes the abundance by a factor **0.4162**.

## Historical comparison

The deterministic sample contains 13 stored transition points, 20 h1-resonance
points, 20 h2-resonance points, 20 DD-boundary points and 27 other/exclusion
points. Selection identities and source-file hashes are saved in
`benchmarks/v2/selection.json`.

| Subset/flag | True | False | Unassessed |
|---|---:|---:|---:|
| Experimental | 23 | 77 | 0 |
| DM | 6 | 91 | 3 |
| Desired tree vacuum global | 74 | 26 | 0 |
| Running boundedness | 64 | 35 | 1 |
| Running unitarity | 85 | 15 | 0 |
| Strict theory subset | 47 | 53 | 0 |
| Baryogenesis candidate | 0 | 17 | 83 |
| GW candidate | 3 | 10 | 87 |

The candidate counts include **78 points not eligible for BSMPT**. Of the 22
eligible rows, baryogenesis is unassessed for five and GW for nine. Positive
evidence survives other incomplete calculations. The three unassessed DM
verdicts are outside the new DD table's coverage. There are 96 fully assessed
and four partially assessed zero-temperature records.

Among rows with stored comparable flags, three DM passes become failures and
four failures become passes; each change includes a changed DD verdict. Eight
HiggsBounds passes, 16 HiggsSignals passes and eight W-mass flags become
failures. The table includes old/new numerical diagnostics and the changed
components. Several corrections were applied together, so these are not
one-change causal ablations.

Twenty-one old `thc=True` records fail boundedness when their stored couplings
are checked. Their old source provenance does not establish why the saved
flags passed. These are reported as historical inconsistencies, without
attributing them to the eigensolver change. The 13 archived transition fixtures
had no stored zero-temperature flags; their new results are first assessments.

### Audit benchmarks and point 223615

- The desired vacuum in the deeper-X benchmark is **not global**: the deeper
  stationary minimum has wx = 674.768 GeV and is lower by
  4.57037e9 GeV^4.
- The running counterexample integrates successfully but violates unitarity
  at **594.659631 GeV**. Numerical integration success is not a theory pass.
- The common-pipeline HiggsSignals SM reference is **152.688433539**; its
  decoupling-limit delta-chi-square is zero.
- Point **223615** remains **baryogenesis=False, GW=True**. The critical EW-entry
  jump is approximately 0.136; the large low-temperature field-space jump is
  approximately 6.59. The competing-minimum crossing near 51.55 GeV is bounded
  by **[51.552258, 51.553137] GeV** in the archived common-temperature refinement.
  Separate unresolved ordering near a branch termination around 80.44 GeV
  remains explicit. This does not erase the resolved crossing or the positive
  transition evidence.

The final provenance check found that the first manto validation used BSMPT
3.2.0 while the local build used 3.2.1. The release uses **3.2.1 on both hosts**,
including its temperature-derivative and sound-speed fixes. The replacement
local source is byte-for-byte identical to the source used for all 22 local
thermal assessments. Manto was rebuilt from that canonical source and the
affected point and integrated pilot were repeated. Earlier 3.2.0 results are
preserved outside the release evidence.

## Executables and launch commands

| Runtime | Laptop path | Manto path |
|---|---|---|
| Python | `/Users/apapaefs/Projects/TwoSingletDM/runtime-v2/venv/bin/python` | `/Users/apapaefs/Projects/runtime-v2/venv/bin/python` |
| Default DM | `/Users/apapaefs/Projects/TwoSingletDM/runtime-v2/micromegas_7.1.4/TRSM/main` | `/Users/apapaefs/Projects/runtime-v2/micromegas_7.1.4/TRSM/main` |
| Comparison DM | `/Users/apapaefs/Projects/TwoSingletDM/runtime-v2/micromegas_6.1.15/TRSM/main` | `/Users/apapaefs/Projects/runtime-v2/micromegas_6.1.15/TRSM/main` |
| CalcTemps | `/Users/apapaefs/Projects/TwoSingletDM/runtime-v2/BSMPT-3.2.1/build/bin/CalcTemps` | `/Users/apapaefs/Projects/runtime-v2/BSMPT-3.2.1/build/bin/CalcTemps` |
| MinimaTracer/PhaseProbe | Same directory as CalcTemps | Same directory as CalcTemps |

On manto, from the deployed repository:

```bash
cd /Users/apapaefs/Projects/TwoSingletDM
/Users/apapaefs/Projects/runtime-v2/venv/bin/python tools/run_next_scan.py \
  --config config/next-scan-v2.json \
  --campaign-dir /Users/apapaefs/Projects/TwoSingletDM/output/next-scan-v2-20260923 \
  --run
```

On the laptop:

```bash
cd /Users/apapaefs/Projects/TwoSingletDM/twosingletDM
/Users/apapaefs/Projects/TwoSingletDM/runtime-v2/venv/bin/python tools/run_next_scan.py \
  --config config/next-scan-v2.json \
  --campaign-dir /Users/apapaefs/Projects/TwoSingletDM/twosingletDM/output/next-scan-v2-20260923 \
  --run
```

Omit `--run` to print the exact underlying campaign command. Add `--pilot` and
choose a fresh directory for four draws across two workers. Individual seeds
can be resumed using the actual main-output path in their output manifest;
keep the same source/runtime receipt and settings. A changed physics fingerprint
requires a new campaign or explicit reevaluation.

The automatic writable CalcHEP directories are isolated per worker. Existing
campaign installations and output directories are preserved. Upstream source
locations for the runtime receipts are `../BSMPT` and `../higgstools` locally,
and `/Users/apapaefs/Projects/runtime-v2/BSMPT-upstream-32ccb11e` and
`/Users/apapaefs/Projects/runtime-v2/higgstools` on manto.

## Interpretation limits retained in the release

The tree vacuum and one-loop RG checks define optional subsets, with no
threshold matching. Equilibrium ordering is among the traced basins; it is
distinct from a completed cosmological transition history. The Z2 freeze-out
flag and fixed-mixing VEV/coupling diagnostics do not solve thermal Boltzmann
evolution. Baryogenesis/GW flags are the agreed qualitative candidate screens.
Unavailable calculations, uncovered DD masses and incomplete thermal histories
remain unavailable. The validation sample uses Tmax = 300 GeV; the supplied next
scan uses Tmax = 1000 GeV.

The [five-family plot index](../validation/v2/release/plots/v2-index.html) shows
these distinctions, including missing-result counts and unresolved brackets.
