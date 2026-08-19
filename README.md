# twosingletDM
TRSM + Dark Matter + ElectroWeak Baryogenesis vs. Higgs Boson Pair production

# Instructions:

## Download MG5_aMC and prepare generated processes

- Copy `loop_sm_twoscalar_generic` into the MG5 `models` directory.
- Launch MG5 and generate the process: 
```
./bin/mg5_aMC
import model loop_sm_twoscalar_generic
generate g g > h h [noborn=QCD]
output gg_hh_twoscalar
launch
```
- Enter and proceed to the next screen, then edit the run card for the desired beam energies.
- Edit `MG5_aMC/input/mg5_configuration.txt`, changing:
```
automatic_html_opening = False
```
- By default `generate_mg5_trsm_xsecs.py` uses
  `../MG5_aMC_v3_5_15` relative to this repository. Override that location
  without editing code by setting `TRSM_MG5_LOCATION`.
- The associated-production interface expects generated process directories
  named `gg_heta0` for `g g > h eta0` and `pp_eta0Z` for
  `p p > eta0 z`. The process-to-directory mapping is in `ProcLocation` in
  `generate_mg5_trsm_xsecs.py`.

To run both associated-production processes during a scan, use:

```bash
python3 generate_trsm_points.py 123 \
  --nrandom 500 \
  --run-mg5
```

With `--run-mg5` and no `--mg5-process` arguments, the defaults are
`gg_heta0` and `pp_eta0Z`. A subset can be selected by repeating the option,
for example `--mg5-process gg_heta0`. By default MadGraph is called only for
fully viable points: `evo`, `thc`, `hb`, `hs`, `ewpo`, `wmass`, and aggregate
`dm` must all be `True`.

To drop only the DM requirement while retaining every evolution, theory, and
experimental constraint, use:

```bash
python3 generate_trsm_points.py 123 \
  --nrandom 500 \
  --run-mg5 \
  --mg5-without-dm
```

These non-DM-viable points are retained in the main scan output even when
`dm=False`. Every requested MG5 column is written for all retained rows, using
`nan` for rows that were not eligible, so the TSV header remains stable.
The optional-mode run tag contains `-noDM` to avoid colliding with a fully
viable MG5 scan using the same date and seed.

The interface updates the generated process's parameter card through MadEvent
`set` commands. It supplies `Meta=M2`, `Miota=M3`, all three physical widths,
`k1`, `k2`, `k3`, and the scalar couplings including `kap133` and `kap233`.
Raw rates are stored as `mg5_xsec_gg_heta0_pb` and
`mg5_xsec_pp_eta0Z_pb`. The scan also stores

```text
mono_higgs_xsec_pb = sigma(gg -> h eta0) * BR(eta0 -> iota0 iota0)
mono_z_xsec_pb     = sigma(pp -> eta0 Z) * BR(eta0 -> iota0 iota0)
```

with `h=H1`, `eta0=H2`, and stable `iota0=H3`.

## Get HiggsTools: https://gitlab.com/higgsbounds/higgstools.git and compile it:
in the HiggsTools directory:
```
mkdir build; cd build; cmake ..; make -j12
```
- You may need to ```make install```, otherwise make sure HiggsTools is in your PYTHONPATH
- You also need to place the HiggssBounds and HiggsSignals datasets in the twosingletDM directories: https://gitlab.com/higgsbounds/hbdataset and https://gitlab.com/higgsbounds/hsdataset.
## Execute:
- To begin the default random scan, execute:
```bash
python3 generate_trsm_points.py SEED --nrandom 500
```
where `SEED` is an integer used as the random-number seed. If `--nrandom` is
omitted, the script defaults to 100 random points.

The default samples `M2` uniformly over
`[m2_min, min(m2_max, m3_max - mhiggs)]`, then samples `M3` over
`[max(m3_min, M2 + mhiggs), m3_max]`. This keeps both masses inside their
configured ranges while enforcing the conditional `M3 >= M2 + mhiggs`
hierarchy. To sample both masses independently over their full configured
ranges, use:

```bash
python3 generate_trsm_points.py 123 \
  --nrandom 500 \
  --independent-m3
```

This fills the rectangular `(M2, M3)` scan range and can therefore include
`M3 < M2` as well as `M3 < M2 + mhiggs`. The mode is mutually exclusive with
the exact and approximate resonant-DM mass modes.

Every random-scan range can be overridden on the command line. Any bound that
is omitted keeps the value in the `define ranges here` block of
`generate_trsm_points.py`. For example, a light-DM rectangular scan can use:

```bash
python3 generate_trsm_points.py 97999 \
  --nrandom 10000 \
  --independent-m3 \
  --m2-min 4 --m2-max 1000 \
  --m3-min 4 --m3-max 65 \
  --scan-k133-k233-log \
  --k133-pow-min -4 --k133-pow-max 3 \
  --k233-pow-min -3 --k233-pow-max 5
```

The complete set of optional overrides is:

```text
--m2-min/--m2-max                 --m3-min/--m3-max
--vs-min/--vs-max                 --k1-min/--k1-max
--lx-min/--lx-max                 --lphix-min/--lphix-max
--lsx-min/--lsx-max               --k133-min/--k133-max
--k233-min/--k233-max             --k133-pow-min/--k133-pow-max
--k233-pow-min/--k233-pow-max
```

The `lphix`, `lsx`, and linear `K133`/`K233` bounds retain the established
sampling convention: they bound the first uniform factor, which is then
multiplied by an independent uniform value in `[-1,1]`. The `*-pow-*` options
are base-10 exponent bounds for the signed logarithmic coupling scan. In
contrast, `--m2`, `--m3`, `--vs`, and the other unsuffixed parameter options
select one explicit point; they are not scan-range settings. Resolved bounds
are saved in scan metadata and the checkpoint fingerprint. A resumed campaign
therefore restores its original ranges, and range flags cannot be changed on
`--resume-from`.

In the `vx=0` branch, the
stable `h3` width and visible branching fractions are recorded as zero,
including below the 20 GeV lower edge of the SM Higgs tables. Whenever
`2*M3 < M1` or `2*M3 < M2`, the generator includes the corresponding
`h1/h2 -> h3 h3` invisible width in the physical total width supplied to
HiggsTools and registers it as `HP.Decay.directInv`. The visible branching
fractions passed to HiggsTools remain the pre-invisible/base values so that
HiggsTools rescales them exactly once.

When `2*M2 < M1`, the generator likewise includes the previously omitted
`H1 -> H2 H2` partial width in the physical `H1` width and registers that
two-BSM-particle decay with HiggsTools.

The 13.6 TeV LO single-scalar cross sections and scalar-cascade metadata are
written as named columns. The exclusive one-invisible topology uses

```text
sigma(parent) * BR(parent -> daughter daughter)
              * 2 * BR(daughter -> H3 H3) * (1 - BR(daughter -> H3 H3))
```

and is stored for both `H2 -> H1 H1` and `H1 -> H2 H2`. The factor of two
counts the two assignments of which identical daughter decays invisibly.

For scalar masses in the scan range below the legacy branching-ratio table,
`4 <= M2 < 20 GeV`, the base SM branching fractions and total width come from
the tracked HiggsTools YR4 `SMHiggs` grid and use the same linear interpolation
as HiggsTools. The established legacy cubic interpolation remains unchanged at
and above 20 GeV.

Every new scan also writes a JSON provenance sidecar beside the tab-separated
data file. For example,

```text
output/trsm_points_13.6-YYYYMMDD-888-False_vxzero.dat
output/trsm_points_13.6-YYYYMMDD-888-False_vxzero.metadata.json
```

The sidecar records the seed, requested point count, stopping and output
selection rules, complete generator invocation, parsed options, mass- and
portal-sampling modes, fixed parameters, and the configured and effective
support of every scanned variable. The distinction matters for the default
mass sampler: although the configured range is `M2 = 4--1000 GeV`, its
effective range is `4--875 GeV`, and the `M3` lower endpoint depends on the
sampled `M2`. With `--independent-m3`, the effective ranges are the full
configured rectangle, `M2 = 4--1000 GeV` and `M3 = 65--1000 GeV`.

### Resume an interrupted random scan

Random scans are checkpointed automatically. The adjacent
`<scan-stem>.checkpoint.json` records the state after the last committed raw
draw, including the Python random-number state, raw draw index, constraint
counters, output sizes, and immutable scan-configuration fingerprint. The
default is to checkpoint every completed draw; use `--checkpoint-every N` to
change that cadence.

Resume by naming the existing main TSV:

```bash
./trsmdm/bin/python generate_trsm_points.py \
  --resume-from output/trsm_points_13.6-20260723-879-False_vxzero.dat \
  --nrandom 10000
```

On resume, `--nrandom` is the **total campaign target**, not the number of
additional points. For a campaign using `--nrandom-count-evo-thc`, the example
therefore stops when the file contains 10,000 evo/thc-passing points. Omitting
`--nrandom` retains the most recently recorded target. The original seed,
sampling ranges and modes, output selection, Higgs/micrOMEGAs conventions, and
BSMPT executable and work-directory settings are loaded automatically from the
metadata sidecar; do not repeat them. They are immutable within one campaign.
The only resume overrides are `--nrandom`, `--checkpoint-every`,
`--print-info`, and `--no-print-info`.

The generator uses an adjacent lock file to reject concurrent writers. A
same-host lock is archived only when its PID is no longer alive. If a crash
leaves output beyond the last checkpoint, the uncommitted tail is backed up and
rolled back before continuation. BSMPT `point_XXXXXX` directories beyond the
committed raw draw are renamed with an `.interrupted-...` suffix rather than
deleted.

Scans created before checkpoint support can be adopted once when both
`--nrandom-count-evo-thc` and `--write-evo-thc-points` were used. Stop the old
generator first and confirm that its PID has exited, then use the same
`--resume-from` command. The generator replays only the inexpensive random
candidate draws from the saved seed and matches every stored scan coordinate
in order. HiggsTools, micrOMEGAs, and BSMPT are not called during this
reconstruction. Adoption aborts before appending anything if the recorded
prefix cannot be reproduced exactly.

The metadata retains the initial target and invocation, adds an append-only
resume history, and reports the current campaign status and counters. A
requested target below the completed count is rejected; an equal target is a
clean no-op marked complete, and a larger target extends a completed campaign.
Starting a fresh command over an existing checkpointed run tag is also rejected
with an instruction to use `--resume-from`.

By default, `--nrandom` is the number of random draws attempted. To instead
keep drawing until `--nrandom` points have passed the RGE evolution (`evo`) and
theory-constraint (`thc`) checks, use:

```bash
python3 generate_trsm_points.py 123 \
  --nrandom 500 \
  --nrandom-count-evo-thc
```

This option only changes the stopping condition. Later constraints such as
HiggsTools, EWPO, W mass, dark matter, and EWPT still determine whether a point
is written or counted as fully viable.

The W-mass constraint uses the tabulated two-state singlet limit in
`datafiles/Tania_MW_SnowmassWhitepaper.dat`. It is applied to
`abs(sin(a12))` and is defined only for `133 <= M2 <= 999` GeV. Points outside
that table range pass the W-mass viability check with a runtime warning because
the constraint is not applicable there; the code does not extrapolate into an
unsupported mass region.

To write those `evo`/`thc`-passing points even when later constraints fail, use
`--write-evo-thc-points`:

```bash
python3 generate_trsm_points.py 123 \
  --nrandom 500 \
  --approximate-resonantDM \
  --delta-res 10 \
  --nrandom-count-evo-thc \
  --write-evo-thc-points
```

This writes points with `evo=True` and `thc=True` to the main output file. It is
different from `--write-all-points`, which writes every evaluated point,
including points that fail `evo` or `thc`.

By default, the vx=0 random scan samples `lphix` and `lsx` directly. To scan
instead over the physical dimensionful couplings `K133` and `K233`, use:

```bash
python3 generate_trsm_points.py 123 \
  --nrandom 500 \
  --scan-k133-k233
```

The `K133` and `K233` scan ranges are set in the `define ranges here` block of
`generate_trsm_points.py`:

```python
K133_min = 1E-6
K133_max = 8.0

K233_min = 1E-4
K233_max = 8.0
```

To scan signed `K133` and `K233` magnitudes logarithmically, use:

```bash
python3 generate_trsm_points.py 123 \
  --nrandom 500 \
  --scan-k133-k233-log
```

The logarithmic scan samples powers set by:

```python
K133_pow_min = -4
K133_pow_max = 3

K233_pow_min = -3
K233_pow_max = 5
```

In either K-scan mode, each sampled `K133` and `K233` point is converted to the
potential-basis inputs expected by the existing BSMPT and micrOMEGAs pipeline:

```text
lphix = 2 * (cos(a12) * K133 + sin(a12) * K233) / 246
lsx   = 2 * (-sin(a12) * K133 + cos(a12) * K233) / vs
```

The canonical `vx=0` convention is

```text
K133 = (lphix * v * cos(a12) - lsx * vs * sin(a12)) / 2
K233 = (lphix * v * sin(a12) + lsx * vs * cos(a12)) / 2
```

where `K_i33` is the coefficient of `h_i h3 h3` in the scalar potential. The
corresponding identical-scalar partial width is
`Gamma(h_i -> h3 h3) = K_i33^2 * beta / (8*pi*M_i)`, and is zero at or below
threshold. Scan files using this normalization carry the convention identifier
`trsm_vxzero_canonical_v1` in both `portal_convention` and
`micromegas_model_convention`.

The scan output includes both the existing `K133` column and a derived `K233`
column so the sampled physical couplings can be inspected directly.

To print the full point table, constraint summary, and dark-matter summary for
each fully evaluated vx=0 scan point, add `--print-info`:

```bash
python3 generate_trsm_points.py 123 \
  --nrandom 500 \
  --write-dm-failed \
  --print-info
```

This is useful together with `--write-dm-failed`: points that pass the non-DM
checks but fail the dark-matter check are still written to the `_dm_failed`
sidecar, and their printed diagnostics include the DM failure reason.
During random scans, `--print-info` also prints a progress counter after each
evaluated point. In `--nrandom-count-evo-thc` mode this includes the number of
random draws, the number of `evo`/`thc`-passing points collected toward the
target, and the number of fully viable points.

- To evaluate one explicit vx=0 point, provide the same point parameters used by
the EWPT helper:
```bash
python3 generate_trsm_points.py 123 \
  --m2 380 --m3 500 --vs 200 --a12 -0.15 \
  --lx 0.10 --lphix 0.050 --lsx 0.15
```
The random scan remains the default; explicit point mode is used only when all
of `--m2`, `--m3`, `--vs`, `--a12`, `--lx`, `--lphix`, and `--lsx` are given.
The current vx=0 random scan samples both positive and negative `a12`.

For resonant dark-matter scans, `--resonantDM1` and `--resonantDM2` set `m3`
from a mass relation instead of sampling or requiring it:

```text
--resonantDM1: m3 = m1/2 = 62.545 GeV
--resonantDM2: m3 = m2/2
```

The flags are mutually exclusive. In explicit point mode, either resonant flag
means `--m3` can be omitted because it is computed from the chosen relation:

```bash
python3 generate_trsm_points.py 123 \
  --m2 380 --vs 200 --a12 -0.15 \
  --lx 0.10 --lphix 0.050 --lsx 0.15 \
  --resonantDM2 \
  --print-info
```

In a random scan, `--resonantDM1` fixes `m3` to the SM-like Higgs-resonant
value for all points, while `--resonantDM2` updates `m3` point-by-point after
each random `m2` draw. In both cases the annihilating dark-matter pair obeys
`2*M3 = M1` or `2*M3 = M2`, respectively.

To scan uniformly near either mass-doubling relation, use
`--approximate-resonantDM` with a mass-window half-width:

```bash
python3 generate_trsm_points.py 123 \
  --nrandom 500 \
  --approximate-resonantDM \
  --delta-res 10
```

For each random point, this mode chooses one of the two approximate branches and
samples uniformly inside the available scan range:

```text
M2 = 2*M3 +/- delta_res
M3 = 2*M2 +/- delta_res
```

If both branches have support in the configured M2/M3 ranges, each is selected
with 50% probability. If only one branch has support, that branch is used for
every point; the scan fails early if neither branch is available. The
approximate mode is mutually exclusive with `--resonantDM1` and
`--resonantDM2`. The `--delta-res` value is the half-width in GeV.

- To run BSMPT EWPT checks only after a generated point passes the existing
viability checks, add `--run-ewpt`:
```bash
python3 generate_trsm_points.py 123 \
  --nrandom 500 \
  --write-dm-failed \
  --run-ewpt \
  --run-ewpt-on-dm-failed \
  --ewpt-require-eq418 \
  --ewpt-thigh 1000 \
  --ewpt-plot-phases \
  --ewpt-workdir ../tests/trsm-ewpt-seed-123
```
`--ewpt-require-eq418` applies the Eq. 4.18 quartic-coupling positivity check
before launching BSMPT, which avoids spending time on EWPT runs that fail this
analytic prefilter. Viable candidate information is printed before this skip, so
discarded points can still be inspected later.

When `--run-ewpt` is enabled and BSMPT returns a first-order transition
strength, the viable-point TSV also includes:

```text
ewpt_ew_true_over_T
ewpt_ew_jump_over_T
```

Both values describe the same selected BSMPT transition. It is chosen from the
available strengths with priority `nucl`, then `perc`, then `compl`, then
`crit`; within one temperature kind the largest finite `ew_true/T` is selected.
If EWPT is not run, Eq. 4.18 skips the run, or no finite strength is available,
the columns are written as `nan`.

`--write-dm-failed` writes points that pass the non-DM checks but fail the
dark-matter check to a separate sidecar file:

```text
output/trsm_points_<run-tag>_dm_failed.dat
```

The usual `trsm_points_<run-tag>.dat` file remains reserved for points that pass
the full viability selection.

The dark-matter check combines the relic-density upper bound, the rescaled
direct-detection limit, and the Fermi-LAT R16 gamma-line indirect limit parsed
from `FermiLAT_line_channel` lines in the micrOMEGAs output. Scan files include
the corresponding diagnostics:

```text
dm_indirect_available
dm_indirect_energy
dm_indirect_flux
dm_indirect_limit
dm_indirect_ratio
dm_indirect_detection_excluded
```

`dm_indirect_ratio` is `dm_indirect_flux / dm_indirect_limit`; values above 1
fail the dark-matter check. For an underabundant candidate with
`xi = Omega / 0.121`, the direct-detection rate is rescaled by `xi` and the
annihilation line flux by `xi^2` (implemented equivalently by dividing the
experimental line-flux limit by `xi^2`).

`--run-ewpt-on-dm-failed` is an exploratory option for otherwise-good points
that fail only the dark-matter check. It runs BSMPT for those points and writes
them to the `_dm_failed` sidecar, including `ewpt_ew_true_over_T` and
`ewpt_ew_jump_over_T` when BSMPT returns a finite strength. Use it together
with `--run-ewpt` if you want BSMPT for both viable and DM-failed points; by
itself it only targets DM-failed points.

### Run EWPT on points from a completed scan without requiring DM

`reprocess_trsm_ewpt.py` runs BSMPT on rows already stored in a completed scan;
it does not regenerate points, rerun HiggsTools or micrOMEGAs, or invoke MG5.
The source TSV is never modified. A row is selected when all stored non-DM
flags pass,

```text
evo & thc & hb & hs & ewpo & wmass
```

while its stored aggregate `dm` result is preserved and reported but does not
gate EWPT. This is appropriate for a source generated with
`--write-evo-thc-points`, which retains every possible non-DM-viable candidate.

For the seed-27999 scan, replace `YYYYMMDD` below with the date in the actual
input filename:

```bash
./trsmdm/bin/python reprocess_trsm_ewpt.py \
  output/trsm_points_13.6-YYYYMMDD-27999-True_vxzero.dat \
  --output output/trsm_points_13.6-YYYYMMDD-27999-True_vxzero_ewpt_non_dm.dat \
  --ewpt-workdir output/ewpt_seed_27999_reprocessed_non_dm \
  --ewpt-thigh 1000 \
  --ewpt-executable /Users/apapaefs/Projects/BSMPT/build/macos-armv8-release/bin/CalcTemps \
  --ewpt-minima-executable /Users/apapaefs/Projects/BSMPT/build/macos-armv8-release/bin/MinimaTracer
```

All input rows and columns, including DM and MG5 results, are copied to the new
TSV. The seven EWPT summary columns are added or updated only for selected rows.
Rows already containing an EWPT attempt are preserved by default, so a scan
originally run with `--run-ewpt` does not repeat its fully viable BSMPT jobs;
add `--rerun-existing-ewpt` only when they should be recalculated as well.
Individual BSMPT failures are recorded as `ewpt_status=failed` and processing
continues.

Progress is committed after every input row to
`<output>.partial`. After an interruption, repeat the command with `--resume`;
the source checksum, destination, work directory, selection and BSMPT settings
must match. The final TSV appears only after every source row is accounted for.
The command also writes an `*.ewpt-reprocess.json` provenance summary and, when
the source has scan metadata, a matching metadata sidecar for the new TSV.
Use `--checkpoint-every N` to change the transaction size and
`--ewpt-require-eq418` to apply the same analytic prefilter available during
generation.

- If the EWPT campaign logs show `ModuleNotFoundError` for packages such as
`scipy`, run with the same Python interpreter used in the working environment,
for example `/Users/apapaefs/.venvs/compphys/bin/python`.

## Plot TRSM Scan Observables

`plot_trsm_observables.py` makes human-editable scatter plots from the TSV files
written by `generate_trsm_points.py`. The default preset is
`ewpt_ew_true_over_T` vs `M2`:

```bash
python3 plot_trsm_observables.py \
  output/trsm_points_13.6-20260529-1234-False_vxzero_dm_failed.dat \
  --output-dir plots \
  --format both
```

Useful one-off customizations:

```bash
python3 plot_trsm_observables.py output/trsm_points_13.6-20260529-1234-False_vxzero_dm_failed.dat \
  --x M2 \
  --y ewpt_ew_true_over_T \
  --color-by dm_omega \
  --size-by dm_dir_det \
  --marker-by dm_relic_excluded \
  --output-stem ewpt_vs_M2_dm_style
```

To add a named plot permanently, edit the `PLOT_PRESETS` dictionary near the top
of `plot_trsm_observables.py`. To inspect available columns:

```bash
python3 plot_trsm_observables.py output/trsm_points_13.6-20260529-1234-False_vxzero_dm_failed.dat --list-columns
```

The plotting script also supports derived observables, defined as Python
functions in the `DERIVED_OBSERVABLES` dictionary near the top of the file. For
example:

```python
DERIVED_OBSERVABLES = {
    "M2_over_M3": lambda row: safe_divide(obs(row, "M2"), obs(row, "M3")),
    "log10_dm_omega": lambda row: safe_log10(obs(row, "dm_omega")),
}
```

Derived observables are added as columns automatically and can be used anywhere
raw columns can be used:

```bash
python3 plot_trsm_observables.py output/trsm_points_13.6-20260529-1234-False_vxzero_dm_failed.dat \
  --x M2_over_M3 \
  --y ewpt_ew_true_over_T \
  --color-by log10_dm_omega \
  --output-stem ewpt_vs_M2_over_M3
```

### Plot the full TRSM constraint suite

For a scan-level overview of the dark-matter and experimental constraints, run
the constraint suite from the repository root:

```bash
/opt/homebrew/bin/python3.11 plot_trsm_constraint_suite.py \
  output/trsm_points_NEW.dat
```

The command above uses the Python installation tested on `manto`. With the
defaults, the suite writes PNG and PDF versions of up to 69 standalone figures
and twelve combined dashboards, together with `constraint_summary.tsv`, under
`plots/trsm_points_NEW_constraints/`. It also writes a self-contained
`index.html` with dashboard and individual-plot previews, links to every
generated PNG/PDF, skipped-plot notices, and the constraint-summary table. Open
that file in a browser to browse the complete suite. When the matching
`<input-stem>.metadata.json` sidecar is present, the index also shows the seed,
generator modes, exact configured/effective variable ranges, fixed parameters,
and full invocation. It shows observed extrema of the retained rows in a
separate column because those can be narrowed by selections. For a legacy scan
without a sidecar, the index labels those extrema explicitly as observed-only
ranges rather than configured bounds. The output location, format, and raster
resolution can be changed explicitly:

```bash
/opt/homebrew/bin/python3.11 plot_trsm_constraint_suite.py \
  output/trsm_points_NEW.dat \
  --output-dir plots/my_new_scan_constraints \
  --format both \
  --dpi 200
```

If a scan or its metadata sidecar was renamed, supply the provenance file
explicitly:

```bash
/opt/homebrew/bin/python3.11 plot_trsm_constraint_suite.py \
  output/trsm_points_RENAMED.dat \
  --scan-metadata output/original_scan.metadata.json
```

The input must be a tab-separated scan file with the standard column names
written by `generate_trsm_points.py`. Selection booleans are parsed strictly as
`True` or `False`. The detailed DM component flags additionally accept the
canonical `nan` written when the microOMEGAs evaluation fails; such relic,
direct-detection, and indirect-detection results are shown as unavailable and
are never treated as passing. Their points still retain the strict aggregate
`dm=False` selection. The suite uses the stored constraint flags rather than
reconstructing the selections. It defines the experimental selection as
`hb & hs & ewpo & wmass`, the theory selection as `evo & thc`, and full
viability as theory, experimental, and dark-matter (`dm`) selections passing
together.

The `23_h2_width_dm_status_m2_m3` map uses the same aggregate DM pass/fail
marker distinction as plot 02, while point color shows the stored complete
$h_2$ width `w2` on a logarithmic scale.

Plots 24--29 reproduce the collaborator diagnostic set as cumulative survivor
overlays for $M_2$--$M_3$, $M_2$--$v_s$, $M_2$--$a_{12}$,
$M_3$--$\lambda_X$, $M_3$--$\lambda_{\Phi X}$, and
$M_3$--$\lambda_{SX}$. The layers are nested in the exact supplied order:
all stored rows, HB, HB and HS, HB and HS and the $W$-mass constraint, and
finally those three constraints together with aggregate DM passing. Distinct
markers and the suite's colorblind-safe palette make the smaller survivor sets
visible above the full scan. This collaborator sequence deliberately does not
apply EWPO and is therefore separate from the suite's `experimental` and
`full_viability` definitions. The fourth dashboard collects all six plots, and
`constraint_summary.tsv` records each cumulative count explicitly.

When the scan contains recorded BSMPT results, plots 30--36, the additional
plot 31b, and a fifth dashboard are added automatically:

- BSMPT run/failed/no-selected-FOPT/weak-FOPT/strong-FOPT status on the
  \(M_2,M_3\) plane;
- the selected \(v_{\rm EW,true}(T_*)/T_*\) and
  \(\Delta v_{\rm EW}(T_*)/T_*\) on the mass plane, each with a
  threshold-centred normalization;
- MinimaTracer global phase-route and electroweak-entry-step maps;
- \(v_{\rm EW,true}(T_*)/T_*\) versus \(M_2\) and \(M_3\); and
- BSMPT evaluation, phase-history, and transition-result counts.

The generator selects \(T_*\) with priority nucleation, then percolation,
completion, and critical temperature. The displayed
\(v_{\rm EW,true}(T_*)/T_*\geq1\) split is the conventional strong-first-order
transition diagnostic; it is not added to the scan's viability selection.
The phase-route map distinguishes direct electroweak entry, paths through an
\(S\)-broken phase, paths through an \(X\)-broken phase, and other multistep
histories. Since an \(X\)-broken phase violates the nominal dark-sector
\(Z_2\), those points are also counted explicitly in `constraint_summary.tsv`.
Rows for which BSMPT was not requested remain visible as “not run”, while
failed evaluations are never interpreted as successful or as having no
first-order transition.

For a file with no recorded BSMPT attempt, these eight figures and their
dashboard are skipped. Cascade/MadGraph figures are independently skipped when
their named scan columns are absent, so legacy inputs remain usable; the HTML
index explains every unavailable plot. A finite `ewpt_ew_true_over_T` is
sufficient for compatibility with older files. Newer files additionally use
`ewpt_ew_jump_over_T`, `ewpt_status`, `ewpt_global_phase_path`,
`ewpt_has_x_broken`, and `ewpt_ew_step_index`. Plot 31b is omitted for legacy
scans without the jump column. Detailed transition temperatures remain in the
per-point `ewpt_result.json` files and are not reconstructed by the scan-table
plot suite.

Plots 37--42 show \(K_{133}\) and \(K_{233}\) separately versus \(M_3\), with
the signed resonance displacement \(M_2-2M_3\) as a symmetric-log color scale
centered on zero. Its dark neutral center makes resonant points visible, while
blue and orange distinguish the two sides. For each coupling the suite writes
variants containing all stored points, only points passing the combined
experimental selection, and only points passing the relic-density constraint.
Relic-density pass requires an available micrOMEGAs relic result with
`dm_relic_excluded=False`; an
unavailable result is never counted as passing. The first new dashboard places
the six variants side by side. Axes and the resonance color normalization are
fixed from all finite stored rows so the three selections are directly
comparable.

Plots 43--46 show the \(M_2\)--\(M_3\) plane colored by the signed \(K_{133}\)
or \(K_{233}\). One variant colors only the combined-experimental survivors and
one colors only aggregate-`dm` survivors for each coupling. Every stored point
is retained as a light-gray reference layer, so the colored selection can be
read against the original scan support. The second new dashboard collects
these four mass-plane maps; each coupling uses a common color normalization
across its experimental and DM variants. Here “experimental” continues to mean
`hb & hs & ewpo & wmass`; the DM maps use the stored aggregate `dm` flag, not
the relic-density component alone.

Plot 47 shows the signed mixing angle `a12` versus `M2`, with the same four-way
DM/experimental categorization used elsewhere in the suite. Plots 40--42
already provide the requested `M3` versus `K233` views for all stored,
experimental-passing, and relic-density-passing points.

Plots 48--51 show the two exclusive one-invisible scalar cascades versus both
`M2` and `M3`, using a logarithmic cross-section axis and the full-viability
selection. Their rates are

```text
H2 -> H1 H1 -> (H1 -> H3 H3) + H1
H1 -> H2 H2 -> (H2 -> H3 H3) + H2
```

with exactly one daughter required to take the invisible branch. Plots 52--55
show the MadGraph mono-Higgs and mono-Z products versus both masses. Those
figures also require full viability. Plots 56--63 repeat the cascade,
mono-Higgs, and mono-Z figures with only the DM requirement removed; all
evolution, theory, and experimental constraints must still pass.
The scalar-cascade and MadGraph results each have a dedicated four-panel
dashboard for each selection. Zero kinematic rates are counted in the
annotation but omitted from the logarithmic y axis; missing or all-`nan` MG
results are reported as unavailable rather than being interpreted as zero.

When the scan stores `k2`, `w2`, and `h2_h3h3_br` and contains at least one
full-viable point with an open `h2 -> h3 h3` decay, plots 64--69 and a twelfth
dashboard are added automatically. They show:

- the `(M2,M3)` plane colored by
  `(sigma_ggF + sigma_VBF) * BR(h2 -> h3 h3)`;
- separate ggF and VBF signal rates versus `M2`, including YR4 production-rate
  uncertainty intervals and a secondary raw-event axis for 3/ab;
- the dominant ggF+VBF rate versus the dark-matter mass `M3`;
- `k2^2` versus `BR(h2 -> h3 h3)`, colored by the signal rate;
- `Gamma2/M2` versus the signal rate, with 1% and 10% width guides; and
- the signal rate versus the direct-detection ratio, colored by the relic
  density ratio.

Only points passing theory, experimental, and aggregate DM constraints are
used in these signal plots. Production rates are calculated as

```text
sigma_P(pp -> h2 -> h3 h3)
  = k2^2 * sigma_P^YR4(M2) * BR(h2 -> h3 h3),  P = ggF, VBF.
```

The tracked table
`datafiles/lhchxswg_yr4_bsm_13p6tev_ggf_vbf.tsv` is extracted from the official
[LHCHXSWG cross-section repository](https://gitlab.cern.ch/LHCHIGGSXS/LHCHXSWG1/crosssections),
workbook `YR4/Higgs_XSBR_YR4.xlsx`, sheet `YR4 BSM 13.6 TeV`, at repository
commit `aad67de39778537fa36fa0692abf66fd43f660a4`. The central ggF and VBF
cross sections cover 10--3000 GeV. The suite interpolates their logarithms
linearly in mass and never extrapolates. It interpolates the separate scale and
PDF+alpha_s percentages linearly, then combines those two components in
quadrature for each displayed uncertainty interval. Cite
[CERN Yellow Report 4](https://arxiv.org/abs/1610.07922) when using these
predictions.

These YR4 BSM cross sections assume the narrow-width approximation and do not
include electroweak corrections. Marker shapes therefore distinguish
`Gamma2/M2 < 1%`, `1% <= Gamma2/M2 < 10%`, and `Gamma2/M2 >= 10%`; the last
category is retained for diagnosis but the factorized NWA rate should not be
used quantitatively there. The plotted rates are inclusive production proxies,
not fiducial search predictions. In particular, an invisible ggF analysis
needs a recoil object such as an ISR jet, while VBF provides an experimentally
direct invisible-Higgs topology. The 3/ab event axis is before acceptance,
triggering, reconstruction, and backgrounds. The HTML index and
`constraint_summary.tsv` record these assumptions, the source provenance, the
available signal-point count, and the width-category counts. Legacy scans
without the required columns remain supported and receive an explicit
signal-section omission notice.

For scans made with `--independent-m3`, the dashed
`M3 = M2 + 125 GeV` line in the mass-plane figures is only a reference to the
default conditional relation; it is not a selection applied to the data in
that mode. The same figures show `M2 = 2*M3` (magenta dash-dot), the physical
threshold for `h2 -> h3 h3`, and the reciprocal mass-hierarchy reference
`M3 = 2*M2` (blue dashed). The latter is not an `h3` decay threshold because
`h3` is the stable dark-matter state. `constraint_summary.tsv` records the
number of points below the first guide and above the second, as well as how
many points have a kinematically open invisible Higgs decay and how many of
those are unmodelled.
Canonical-v1 files should have zero unmodelled points; legacy files without the
new metadata remain supported and are classified conservatively.

For informative pass/fail comparisons, generate the input scan with
`--write-evo-thc-points` to retain all points passing `evo` and `thc`, or use
`--write-all-points` to retain every evaluated point. The generator's default
output contains only fully viable points, so it provides little or no
separation between the constraint categories.

The suite tests can be run independently with:

```bash
/opt/homebrew/bin/python3.11 -m unittest test_plot_trsm_constraint_suite.py
```

### Re-evaluate a legacy scan with canonical invisible widths and DM normalization

Legacy scan files are never modified in place. To recompute the physical Higgs
widths, HiggsBounds/HiggsSignals details, and all micrOMEGAs result fields with
the canonical-v1 conventions, run:

```bash
trsmdm/bin/python reevaluate_trsm_dm_higgs.py \
  output/trsm_points_OLD.dat \
  --output output/trsm_points_OLD_canonical_v1.dat \
  --micromegas-main ../micromegas_6.1.15/TRSM/main \
  --checkpoint-every 25
```

The calculation is sequential because each point invokes the shared
micrOMEGAs executable. Progress is checkpointed to an adjacent `.partial`
file. If a run is interrupted, repeat the command with `--resume`; the script
checks that the partial file belongs to the same input before continuing.
Existing theory, EWPO, W-mass, EWPT, MG5, and unknown columns are preserved.
The output path must differ from the input and an existing completed output is
not overwritten.

After completion, regenerate the constraint suite from the versioned file:

```bash
/opt/homebrew/bin/python3.11 plot_trsm_constraint_suite.py \
  output/trsm_points_OLD_canonical_v1.dat
```

Historical DM columns in legacy scans use the old CalcHEP normalization and
must not be corrected by numerical rescaling; they need a real re-evaluation.
The separate electroweak-input difference between the generated CalcHEP value
of `v` and the Python value `v=246 GeV` is deliberately left isolated and is
not folded into the portal-coupling normalization.

## Run TRSM EWPT checks with BSMPT

The helper script `test_trsm_ewpt.py` writes a one-point `TRSM_Input.tsv`, runs
BSMPT `MinimaTracer` first, reconstructs the global minimum branch, and then runs
`CalcTemps`. The default input basis is

```text
m1 m2 m3 vs a12 lx lphix lsx
```

where `m1` is fixed to 125.09 GeV by the script, `m2` and `m3` are scalar
masses, `vs` is the singlet vev, `a12` is the doublet-singlet mixing angle, and
`X` has zero zero-temperature vev.

The examples below assume they are run from this `twosingletDM` repository
directory. From the parent `TwoSingletDM` project directory, prefix the script
path with `twosingletDM/` and use `tests/...` instead of `../tests/...`.

The script defaults to the local BSMPT build used in this checkout:

```text
/Users/apapaefs/Projects/TwoSingletDM/BSMPT/build/macos-armv8-release/bin/CalcTemps
```

Override the binaries if needed:

```bash
python3 test_trsm_ewpt.py \
  --executable /path/to/CalcTemps \
  --minima-executable /path/to/MinimaTracer \
  --m2 160 --m3 420 --vs 70 --a12 0.20 \
  --lx 0.12 --lphix 0.03 --lsx 0.10
```

Example run with MinimaTracer phase plots:

```bash
python3 test_trsm_ewpt.py \
  --m2 160 --m3 420 --vs 70 --a12 0.20 \
  --lx 0.12 --lphix 0.03 --lsx 0.10 \
  --thigh 1000 \
  --plot-phases \
  --plot-output phases \
  --workdir ../tests/trsm-ewpt-example
```

Plotting requires `matplotlib`; without `--plot-phases`, the text and JSON
summaries do not need it.

This prints the `CalcTemps` output path, the `MinimaTracer` output path, optional
PNG/PDF phase plots, the compressed cooling history, and transition strengths,
for example:

```text
global_phase_path: SYM -> SINGLET_S -> EW
status_nucl_0: success
fopt_strength_nucl_0: T=..., ew_jump/T=..., ew_true/T=..., field_jump/T=...
```

For a two-step singlet-assisted EWPT search, look for a cooling path like

```text
SYM -> SINGLET_S -> EW
```

and then require the EW step to complete, e.g. `status_nucl_0: success`. The
critical-temperature result alone is not enough if `status_bounce_sol_0` fails
or `T_nucl_0` is `nan`. As a rough baryogenesis diagnostic, inspect
`ew_jump/T` or `ew_true/T`; values near or above 1 are the interesting regime.

A light-singlet, paper-inspired example scan point is:

```bash
python3 test_trsm_ewpt.py \
  --m2 5 --m3 600 --vs 500 --a12 -0.25 \
  --lx 0.25 --lphix 0.001 --lsx 0.001 \
  --thigh 1000 \
  --plot-phases \
  --plot-output phases \
  --workdir ../tests/trsm-ewpt-1911-m2-5-vs-500-a12m025-xdecoupled
```

Small scans can be done directly from the shell:

```bash
for M2 in 5 8 10; do
  for VS in 220 500 800; do
    python3 test_trsm_ewpt.py \
      --m2 "$M2" --m3 600 --vs "$VS" --a12 -0.20 \
      --lx 0.25 --lphix 0.001 --lsx 0.001 \
      --thigh 1000 \
      --plot-phases \
      --plot-output phases \
      --workdir "../tests/trsm-ewpt-m2-${M2}-vs-${VS}-a12m020"
  done
done
```

Use `--json` to capture the first `CalcTemps` row plus the nested MinimaTracer
phase summary:

```bash
python3 test_trsm_ewpt.py \
  --m2 10 --m3 600 --vs 500 --a12 -0.20 \
  --lx 0.25 --lphix 0.001 --lsx 0.001 \
  --thigh 1000 \
  --json
```

Run the focused test suite after editing the EWPT runner:

```bash
python3 test_trsm_ewpt_runner.py
```

## Run Multi-Seed TRSM Campaigns

`run_trsm_seed_campaign.py` launches `generate_trsm_points.py` over a contiguous
seed range, keeps a live per-seed log, aggregates viable point files, and ranks
the best EWPT points by `ew_jump/T`.

Example from this `twosingletDM` directory:

```bash
python3 run_trsm_seed_campaign.py \
  --seed-start 1 \
  --nseeds 20 \
  --nrandom 500 \
  --jobs 4 \
  --heartbeat-seconds 30 \
  --campaign-dir ../tests/trsm-campaign-001 \
  --run-cwd /Users/apapaefs/Projects/TwoSingletDM/twosingletDM \
  --python-executable /Users/apapaefs/.venvs/compphys/bin/python \
  --write-dm-failed \
  --run-ewpt \
  --run-ewpt-on-dm-failed \
  --ewpt-require-eq418 \
  --ewpt-thigh 1000
```

Example from the parent `TwoSingletDM` directory:

```bash
python3 twosingletDM/run_trsm_seed_campaign.py \
  --seed-start 1 \
  --nseeds 20 \
  --nrandom 500 \
  --jobs 4 \
  --heartbeat-seconds 30 \
  --campaign-dir tests/trsm-campaign-001 \
  --run-cwd /Users/apapaefs/Projects/TwoSingletDM \
  --python-executable /Users/apapaefs/.venvs/compphys/bin/python \
  --write-dm-failed \
  --run-ewpt \
  --run-ewpt-on-dm-failed \
  --ewpt-require-eq418 \
  --ewpt-thigh 1000
```

While running, the launcher prints startup, launch, heartbeat, and completion
lines such as:

```text
Starting TRSM seed campaign: seeds=1..20 nseeds=20 jobs=4 logs=...
[launch] seed=1
[running] active=4 seeds=1,2,3,4
[3/20 done] seed=17 viable=4 ewpt_runs=2 best_ew_jump/T=0.31 elapsed=...
```

For detailed live output from one seed, tail its log in another terminal:

```bash
tail -f ../tests/trsm-campaign-001/logs/seed_1.log
```

Campaign outputs are written to:

```text
campaign_summary.tsv
campaign_summary.json
combined_points.tsv
best_points.tsv
logs/seed_<seed>.log
ewpt/seed_<seed>/point_*/ewpt_summary.txt
```

`best_points.tsv` uses the largest available first-order EWPT strength with the
priority `nucl`, then `perc`, then `compl`, then `crit`.

Run the focused tests for the scan and campaign helpers with:

```bash
python3 test_generate_trsm_points_runner.py
python3 test_run_trsm_seed_campaign.py
python3 test_trsm_ewpt_runner.py
python3 test_reprocess_trsm_ewpt.py
python3 -m py_compile generate_trsm_points.py run_trsm_seed_campaign.py reprocess_trsm_ewpt.py
```


## Install micrOMEGAs 6.1.15 and create TRSM:

- ```cd DM``` and untar `micromegas_6.1.15.tar`: ```tar xvzf micromegas_6.1.15.tar```

- ```cd micromegas_6.1.15/``` then follow the READ.ME instructions I and V:
```[g]make``` and ```./newProject TRSM```

- copy the desired model files from DM/models into micromegas. I use h4GOn:
```cp ../models/h4GOn/* TRSM/work/models/```

- Apply the tracked `calcSpectrum` finite-value guard. This makes micrOMEGAs
return an error instead of looping in channel sorting when an earlier failed
relic-density calculation leaves non-finite annihilation weights:
```bash
patch -p1 < ../patches/micromegas-6.1.15-calcspectrum-finite-guard.patch
```

- The checked-in `h4GOn` and `h4GOff` models implement the canonical-v1
normalization: the potential contains `LHX/2`, `LSX/2`, and `LX/4`, while the
input card remains one-to-one (`LHX=lPhiX`, `LSX=lSX`, and `LX=lX`). Do not
compensate by scaling card values. If replacing an older installed model, back
up `TRSM/work/models` and the executable before copying the corrected files.

- In `TRSM/` there is a `main.c` file and a Makefile. Rebuild after installing
the model with:
```bash
make clean
make main=main.c
```
The setup is complete and and micromegas is used by codes in DM/example.

- To test manually, copy the example data point:
```cp DM/data.par micromegas_6.1.15/TRSM/```
```./main data.par```
compare the output with DM/example_test.out

## run micrOMEGAs with a steering code:
- check out DM/example_steer/README.md and create another directory for the study you'd like to perform.
