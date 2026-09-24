# Standalone DM and Planck CMB validation

This example uses the **same v2 cards, output parser and constraint assessment as
production** (`test_trsm_DM.py` and `trsm_cmb.py`). It runs independently of
HiggsTools and BSMPT. Use Python **3.10+**; only plots require matplotlib.

Defaults: corrected **micrOMEGAs 7.1.4**, **Planck CMB enabled**, virtual W/Z
channels enabled, the observed **LZ WS2024** SI table, relic upper bound **0.121**,
and abundance reference **0.12**. A rebuilt v2 driver is required; its capability
query rejects older executables. See the [runtime setup](../../docs/constraints-v2.md)
and [CMB prescription](../planck-cmb.md).

## 1. Reproduce the four benchmarks

Run from the **repository root**. With the installed `runtime-v2` environment:

```bash
export TRSM_RUNTIME_ROOT="$(cd ../runtime-v2 && pwd)"
TRSM_PYTHON="$TRSM_RUNTIME_ROOT/venv/bin/python"

"$TRSM_PYTHON" DM/steer_example/validate_cmb.py \
  --output-dir /tmp/trsm-cmb-check-7 --plot

"$TRSM_PYTHON" DM/steer_example/validate_cmb.py \
  --micromegas-version 6 --output-dir /tmp/trsm-cmb-check-6 --plot
```

Choose a **new output directory** each time. On manto, first
`cd /Users/apapaefs/Projects/TwoSingletDM`; the same commands then apply.
On another installation, use your Python executable and add
`--micromegas-main /absolute/path/to/corrected/TRSM/main`. The version option
selects the reference and default executable; a custom executable's identity is
recorded by path, SHA256 and capabilities.

These are four fixed points per backend, not a scan. `validate_cmb.py` explicitly
enables CMB for **both** versions. It checks solver success, the loop hook in
relic and indirect calculations, the input pole mass, virtual W/Z settings,
raw-output parsing, independent abundance/rescaling arithmetic, and the saved
reference results. It exits nonzero on failure. Cross-host numerical tolerance
is relative `2e-6`, absolute `1e-30`; versions are compared to **their own** references.

Validation on 2026-09-24: all **86 checks per backend** pass locally and on
manto; all saved point diagnostics agree exactly between hosts. The example
regressions and the production DM/CMB tests also pass.

The full-precision [point list](benchmarks/cmb-points.dat),
[reference values and provenance](benchmarks/cmb-reference-v2.json), and native
logs in `benchmarks/raw-7.1.4/` and `benchmarks/raw-6.1.15/` are committed.
The following rounded values are from the corrected 7.1.4 driver, with virtual W/Z on:

| Index | Case | MX [GeV] | Omega h² | Raw CMB ratio | Rescaled CMB ratio | CMB |
|---|---|---:|---:|---:|---:|---|
| 1 | Heavy | 500 | 0.266204 | 0.0105482 | 0.0105482 | pass |
| 2 | Below h1 pole | 62 | 0.00183464 | 0.0285223 | 6.66691e-6 | pass |
| 3 | At h1 pole | 62.545 | 0.00309248 | 9078.91 | 6.02954 | excluded |
| 4 | Light, larger portals | 10 | 0.119135 | 1.68753 | 1.66330 | excluded |

Only point 2 passes the aggregate DM assessment. These fixtures are constraint
benchmarks, not points passing every model constraint. The two backend versions
need not agree: for example, the 10 GeV point gives Omega h² = 0.129264 and a
rescaled CMB ratio of 1.56536 with 6.1.15. The reference records retain that difference.

Inspect:

- `cmb-comparison.tsv`: raw and rescaled CMB ratios alongside relic and freeze-out values.
- `validation.json`: every check and its verdict.
- `cmb_ratios.png`: raw and assessed ratios, exclusion line and unavailable count.
- `OUT_mO_<index>`: untouched native stdout, including `TRSM_PlanckCMB_v1`.
- `ERR_mO_<index>`: native stderr.
- `MO_inp<index>.dat`: actual full-precision card, including the shared SM inputs.

This validates the implementation and reproduces native benchmarks; it is not
an independent Boltzmann solver or Planck likelihood calculation. The production
low-velocity, s-wave approximation is retained.

## 2. Check the arithmetic without a micrOMEGAs installation

From the repository root, with Python 3.10+:

```bash
python3 DM/steer_example/run_scan.py \
  --input DM/steer_example/benchmarks/cmb-points.dat \
  --raw-output-dir DM/steer_example/benchmarks/raw-7.1.4 \
  --output-dir /tmp/trsm-cmb-replay

python3 DM/steer_example/validate_cmb.py \
  --output-dir /tmp/trsm-cmb-replay --check-only
```

For the v6 logs, select `raw-6.1.15` and pass `--micromegas-version 6 --planck-cmb`
to `run_scan.py`, then `--micromegas-version 6` to the checker.
Replay validates parsing and comparisons; it does not recalculate the native rates.

The driver returns a full-abundance ratio to
`p_ann < 3.2e-28 cm^3 s^-1 GeV^-1` (95% CL). The wrapper computes

```text
xi = min(1, Omega_h2 / 0.12)
R_CMB = ratio_raw * xi**2
excluded = R_CMB > 1
```

Equality at one passes. `--no-rescale` sets `xi=1` for CMB and removes abundance
rescaling from DD and gamma lines too. Missing/failed CMB results have a **null**
exclusion verdict, never an assumed pass. `--no-planck-cmb` explicitly disables
this constraint. The standalone v6 runner defaults to CMB disabled, matching
production; use `--planck-cmb` to enable it.

## 3. Run or replay one point

```bash
"$TRSM_PYTHON" DM/steer_example/run_single_point.py \
  --index 1 --lx 0.1 --lhx 0.001 --lsx 0.01 --mx 62.545 \
  --vevs 200 --sint -0.14943813247359922 --mh2 380 \
  --output-dir /tmp/trsm-dm-point
```

Alternatively, supply `--card /path/to/MO_inp1.dat`. For another filename also
supply `--index`. Cards use `SinT = sin(a12)` with `a12` in the principal interval
`[-pi/2, pi/2]`. Nonstandard SM overrides and unknown card keys are rejected
rather than silently discarded. The seven legacy scalar inputs remain accepted;
the common SM inputs are added to the actual card.

To reassess that exact raw output with full abundance:

```bash
python3 DM/steer_example/run_single_point.py \
  --card /tmp/trsm-dm-point/MO_inp1.dat \
  --micromegas-output /tmp/trsm-dm-point/OUT_mO_1 --no-rescale \
  --output-dir /tmp/trsm-dm-point-full-abundance
```

`result_<index>.json` and `.tsv` contain the verdicts, reasons, raw/rescaled CMB
ratios, actual masses and widths, Xf, Tf, and loop-hook counts. JSON uses `null`;
TSV uses `nan` for unavailable values. Native errors preserve logs and return a
nonzero shell status. A successfully calculated but excluded point returns zero.
Old logs lacking solver/CMB evidence remain unassessed.

## 4. Run your own point list or generated cards

The whitespace input format is
`index LX LHX LSX MX vevs SinT Mh2`; blank lines and `#` comments are accepted.
Indices must be unique nonnegative integers. Parameters are written with 17
significant digits. For example:

```bash
python3 DM/steer_example/generate_oks.py \
  --lhx-start 0.01 --lhx-end 0.03 --lhx-step 0.01 \
  --output /tmp/my-dm-points.dat

"$TRSM_PYTHON" DM/steer_example/run_scan.py \
  --input /tmp/my-dm-points.dat --output-dir /tmp/my-dm-results
```

The generator also supports per-parameter `--*-log`, `--*-num-points`,
`--equal-couplings` (`LSX=LHX`) and `--resonant-mx` (`MX=Mh2/2`).
For the light Higgs pole, explicitly use `--mx-start 62.545`.

The existing card/shell workflow remains available, from any working directory:

```bash
python3 DM/steer_example/source/write_mo.py \
  --input /tmp/my-dm-points.dat --output-dir /tmp/my-dm-cards

PYTHON="$TRSM_PYTHON" CARD_DIR=/tmp/my-dm-cards OUTPUT_DIR=/tmp/my-dm-results-cards \
  bash DM/steer_example/run/MOrun.sh
```

`MICROMEGAS_MAIN` overrides the executable for this shell entry point. Extra
arguments such as `--micromegas-version 6 --planck-cmb` are forwarded. Each runner
has its own temporary writable CalcHEP directory. Separate jobs must also use
separate output directories. The optional Condor template needs site-specific
Python/runtime paths and a shared filesystem; no job is submitted automatically.

Batch outputs include all points, including exclusions and failures:

| File | Meaning |
|---|---|
| `results.json`, `results.tsv` | Authoritative headed records with every diagnostic and nullable verdict |
| `metadata.json` | Physics profile, source/executable hashes, inputs, SI table checksum and CMB settings |
| `counts.json` | Counts of assessed passes, exclusions and unavailable results |
| `oks.dat` | Input rows corresponding to the processed output |
| `allall.dat`, `dmexcl.dat`, `dmunassessed.dat` | Aggregate DM pass, exclusion, unassessed subsets |
| `cmbpass.dat`, `cmbexcl.dat`, `cmbunassessed.dat` | Separate CMB subsets; disabled CMB appears in none |
| `cmb_caughtit.dat` | CMB exclusions passing the other DM screens |
| `all_dirpass.dat` | Relic and DD passes, irrespective of gamma/CMB results |

`dm_passed` uses the production aggregation policy. A known exclusion can give
`False` even when another constraint is unavailable; its individual status is
still recorded. Gamma-line availability/coverage remains a separate diagnostic.
A solver failure gives an unassessed aggregate result.

Legacy `scan_results.dat` keeps its original **11** columns (first eight inputs,
MDM, Omega, SI in pb). Only successful solver results enter it. `DM_data_<index>`
and subset `.dat` files keep their **18** columns:

```text
index LX LHX LSX MX vevs SinT Mh2 MDM Omega DirDet DirDetLimit DirDetBaseLimit
IndirAvailable IndirEnergy IndirFlux IndirLimit IndirRatio
```

CMB columns are in the headed JSON/TSV outputs. Legacy positional tables alone
cannot establish assessment status. `source/mO_excluder.py` now accepts the same
arguments as `run_single_point.py`; its old positional numbers omitted CMB and
solver evidence and are rejected with a migration message. Replaying full logs
is the replacement. Historical SI fits require explicit `--limit-model
legacy-output` or `--limit-model lz2025-source`; `--limit-table` selects a normalized
alternative table. No SI table extrapolation is performed.

Output directories are protected against overwriting previous results. This
small example has no campaign resume: saved raw logs can be replayed into a new
directory. The production campaign tools provide fingerprinted resume.
`concat_dat.py` is a legacy numeric-table utility; it does not merge v2 metadata
or authoritative JSON/TSV records.

## 5. Plots and regression tests

```bash
"$TRSM_PYTHON" DM/steer_example/plot/plot_cmb.py /tmp/my-dm-results
"$TRSM_PYTHON" DM/steer_example/plot/steer_plots.py --outdir /tmp/my-dm-results
python3 -m unittest discover -s DM/steer_example/tests -v
```

The former limit-plot entry points remain available and now read the shared
limits (matplotlib replaces gnuplot):

```bash
MO_EXCLUDER_OUTPUT_DIR=/tmp/trsm-limit-plots "$TRSM_PYTHON" \
  DM/steer_example/source/mO_excluder.py --plot-dirdet-limits
MO_EXCLUDER_OUTPUT_DIR=/tmp/trsm-limit-plots "$TRSM_PYTHON" \
  DM/steer_example/source/mO_excluder.py --plot-indirect-limits
```

The plotting driver reads the output's own `oks.dat` and includes the CMB figure.
Existing relic/DD/gamma plots retain their explicitly named subset meanings;
use the CMB plot and headed verdicts for the CMB check. The extra
`0.119 <= Omega <= 0.121` plotting band is not an additional acceptance cut.
