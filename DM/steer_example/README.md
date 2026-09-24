# Standalone micrOMEGAs scans and plots

This directory still supports the general workflow: **generate points → run
micrOMEGAs → retain the full output → apply relic-density, direct-detection and
gamma-line cuts → plot the successive selections**. Planck CMB is an additional
constraint; the four CMB benchmarks are optional tests, not the scan input.
HiggsTools and BSMPT are not needed for this workflow.

## Install once

Use a **full repository clone** and Python **3.10 or newer**. The shared module
`test_trsm_DM.py` is a tracked file at the repository root, not a pip package.
The DM runner itself uses the Python standard library. Plots need NumPy and
matplotlib:

```bash
# From the repository root
python3 -m venv .venv-dm
source .venv-dm/bin/activate
python -m pip install numpy matplotlib
```

Install micrOMEGAs following the complete [DM installation instructions](../README.md).
In particular, supply the extracted release directory to the setup script:

```bash
# From the repository root, after extracting the release
sh DM/setup_micromegas.sh ../runtime-v2/micromegas_7.1.4
python DM/steer_example/check_installation.py
```

`setup_micromegas.sh` builds an executable; it is not an environment activation
script. Sourcing it is tolerated and will no longer close the calling shell,
but without the required argument it prints usage and returns status 2.

For an installation elsewhere, use `--micromegas-main /absolute/path/TRSM/main`
or export `MICROMEGAS_MAIN` with that path. CLI takes precedence. Alternatively
set `TRSM_RUNTIME_ROOT` to the parent of `micromegas_7.1.4`.
`--micromegas-version 6` selects a rebuilt 6.1.15 backend. A v2 run requires the
model, loop helper and driver installed together; an older executable cannot
be made v2 simply by changing its path.

## Generate, run and plot your own scan

The following commands run **inside `DM/steer_example`**, with the Python
environment still active:

```bash
cd DM/steer_example
python check_installation.py

python generate_oks.py --lx-start 0.1 --lhx-start 0.01 --lsx-start 0.02 \
  --mx-start 50 --mx-end 60 --mx-step 5 --vevs-start 300 \
  --sint-start 0.1 --mh2-start 200 --output run/my-points.dat

python run_scan.py --input run/my-points.dat --output-dir output/my-scan
python plot/steer_plots.py --outdir output/my-scan
```

Choose a new output directory for another run. `generate_oks.py` supports all
seven scalar inputs, linear or logarithmic ranges, `--equal-couplings`
(`LSX=LHX`), and `--resonant-mx` (`MX=Mh2/2`). See its `--help` for range options.
An input list has eight whitespace-separated columns:

```text
index LX LHX LSX MX vevs SinT Mh2
```

To run the relic/DD/gamma workflow without CMB, add `--no-planck-cmb` to the
`run_scan.py` command. Version 7 defaults to CMB enabled; version 6 defaults to
it disabled and accepts `--planck-cmb` explicitly. `--no-rescale` disables
abundance rescaling for DD, gamma lines and CMB. `--relic-upper-limit` remains
available; its default is 0.121 and the v2 abundance reference remains 0.12.
The current SI limit is the observed LZ WS2024 table. These physics defaults
differ from historical releases; see [the model change audit](../models/changes-v2.md)
and [the constraint prescription](../../docs/constraints-v2.md).

For a two-variable scan, `steer_plots.py` selects axes automatically or accepts
`--xvar MX --yvar LHX`. It preserves the older relic/DD/gamma plots and adds
`constraint_cutflow.png`: all evaluated points, relic passes, then successive
DD, gamma-line and enabled CMB selections. Unassessed points are counted
separately. With one varying input it plots that input against Omega. The new
cutflow and CMB plots go in the scan output; historical plots retain their
`plot/outplots/` destination. Each plot script also accepts `--output`.
The narrow `0.119 <= Omega <= 0.121` plot is a labelled display band, not an
additional acceptance requirement.

## Keep the original card/shell workflow

These existing entry points and default directories still work:

```bash
# Inside DM/steer_example; run/oks.dat is the original input location
cd source
python write_mo.py
cd ../run
./MOrun.sh
cd ..
python plot/steer_plots.py
```

`write_mo.py` also accepts explicit paths. It can regenerate the same card
directory as before. Changed and stale cards are retained in a
`.previous-cards-*` subdirectory, so `MOrun.sh` sees exactly the current input
list. For a separate named run:

```bash
python source/write_mo.py --input run/my-points.dat --output-dir run/my-cards
CARD_DIR="$PWD/run/my-cards" OUTPUT_DIR="$PWD/output/my-card-scan" \
  PYTHON="$(command -v python)" bash run/MOrun.sh --no-planck-cmb
```

`MICROMEGAS_MAIN`, `CARD_DIR`, `OUTPUT_DIR`, `PYTHON`, and `LOG_FILE` are supported
by `MOrun.sh`. Additional runner options are forwarded. No scan or batch script
submits a cluster job. Every native runner uses its own writable CalcHEP cache.

Batch outputs retain all evaluated points, including exclusions and failures:

| Output | Contents |
|---|---|
| `MOrun.log` | Runner terminal log |
| `OUT_mO_<index>`, `OUT_mO/OUT_mO_<index>` | Complete native stdout, in new and historical locations |
| `ERR_mO_<index>` | Native stderr |
| `MO_inp<index>.dat` | Actual full-precision card |
| `results.json`, `results.tsv`, `result_<index>.json/.tsv` | All diagnostics and individual verdicts |
| `metadata.json`, `oks.dat` | Effective physics configuration, provenance and processed inputs |
| `scan_results.dat` | Original 11 columns, for successful solver results |
| `DM_data_<index>`, `DM_data/DM_data_<index>` | Original 18-column point summaries |
| `relic_pass.dat`, `omexcl.dat`, `relic_strict.dat` | Relic subsets and narrow display band |
| `luxpass.dat`, `luxexcl.dat`, `all_dirpass.dat`, `omgpass_dirfail.dat` | DD and cumulative relic/DD subsets |
| `indirpass.dat`, `indirexcl.dat`, `*_caughtit.dat` | Gamma-line and complementary exclusion subsets |
| `all_indirpass.dat` | Cumulative relic/DD/gamma passes before applying CMB |
| `allall.dat`, `dmexcl.dat`, `dmunassessed.dat` | Aggregate pass, excluded and unassessed subsets |
| `cmbpass.dat`, `cmbexcl.dat`, `cmbunassessed.dat` | Separate CMB subsets |
| `counts.json`, `cutflow.json` | Individual and cumulative selection counts |

The 11-column layout is `index LX LHX LSX MX vevs SinT Mh2 MDM Omega DirDet`.
The 18-column layout appends `DirDetLimit DirDetBaseLimit IndirAvailable
IndirEnergy IndirFlux IndirLimit IndirRatio`. CMB and status fields live in the
headed JSON/TSV products, so positional layouts stay unchanged. JSON uses `null`
and TSV uses `nan` for unavailable values. Failed native calculations return
nonzero; a successfully calculated but excluded point returns zero.
Gamma-line coverage remains a separate diagnostic, following production policy.

## One point, saved logs, and legacy positional calls

```bash
# Inside DM/steer_example
python run_single_point.py --index 1 --lx 0.1 --lhx 0.01 --lsx 0.02 \
  --mx 50 --vevs 300 --sint 0.1 --mh2 200 --output-dir output/one-point

python run_single_point.py --card output/one-point/MO_inp1.dat \
  --micromegas-output output/one-point/OUT_mO_1 --output-dir output/replayed-point

python run_scan.py --input run/my-points.dat --raw-output-dir output/my-scan \
  --output-dir output/replayed-scan
```

Replay accepts flat logs and the original nested `OUT_mO/` directory. It needs
no compiled backend. Old logs without solver-status evidence remain explicitly
unassessed in the v2 products.

The original positional excluder is supported again:

```bash
MO_EXCLUDER_OKS_FILE="$PWD/run/my-points.dat" \
MO_EXCLUDER_OUTPUT_DIR="$PWD/output/legacy-assessment" \
  python source/mO_excluder.py 1 50 0.06 1e-15 50 1e-20
```

Its arguments are `index m_DM Omega DirDet [E_gamma Phi_R16]...`. This historical
interface uses the original LZ2025 fit and abundance reference 0.121, writes the
original subset files, and appends one row per call. It cannot assess solver
success or CMB because those values are absent from the interface. It reports
that limitation and keeps its output separate from v2 results. For current
constraints use complete raw-log replay instead. The limit-plot commands remain:

```bash
python source/mO_excluder.py --plot-dirdet-limits
python source/mO_excluder.py --plot-indirect-limits
```

They plot the current shared limits. `MO_EXCLUDER_OUTPUT_DIR` selects their output.

## Optional benchmarks and tests

The benchmark files exercise the same general runner. They do not restrict
which points or constraints you can evaluate:

```bash
python validate_cmb.py --output-dir output/cmb-check --plot
python -m unittest discover -s tests -v
```

To test parsing without an installation:

```bash
python run_scan.py --input benchmarks/cmb-points.dat \
  --raw-output-dir benchmarks/raw-7.1.4 --output-dir output/benchmark-replay
python validate_cmb.py --output-dir output/benchmark-replay --check-only
```

`validate_cmb.py` checks four points against the saved reference for the selected
backend, with relative tolerance `2e-6` and absolute tolerance `1e-30`.
[Reference data](benchmarks/cmb-reference-v2.json) and [native logs](benchmarks/)
are retained. The [Planck CMB note](../planck-cmb.md) explains the s-wave
approximation and the rescaling `min(1, Omega/0.12)^2`.

## Troubleshooting

- **`ModuleNotFoundError: test_trsm_DM`**: no pip installation of that module is
  needed. From `DM/steer_example`, run `python check_installation.py --no-native`.
  It prints the actual interpreter, repository and module paths. Verify that
  `../../test_trsm_DM.py` exists in the clone. The examples resolve imports from
  their own location and do not require `PYTHONPATH` or a particular current
  directory. `TRSM_REPO_ROOT=/absolute/path/to/twosingletDM` can explicitly select
  the checkout. If a complete clone still fails, retain the full traceback,
  exact command, and this diagnostic output; the current regression tests run
  these entry points with `PYTHONPATH` removed.
- **Setup usage/status 2**: give the extracted micrOMEGAs directory as an argument.
  Setup now runs in a subshell, so it cannot call `exit` on your interactive shell
  or leave `set -e`, `set -u`, or a changed directory behind.
- **Missing backend / capability error**: use the [installation steps](../README.md)
  and `check_installation.py --micromegas-main /path/to/TRSM/main`.
- **Existing output**: choose a fresh output directory. This prevents merging
  different physics settings or duplicating results. Saved raw output remains
  replayable. `concat_dat.py` still combines legacy numeric tables; it does not
  merge v2 provenance or status records.
