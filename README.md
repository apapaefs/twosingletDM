# twosingletDM
TRSM + Dark Matter + ElectroWeak Baryogenesis vs. Higgs Boson Pair production

## Current constraint profile

Use [constraint pipeline v2](docs/constraints-v2.md) and the
[next-scan readiness report](docs/next-scan-readiness-v2.md) for the active
`vx=0` scan. The versioned configuration is
[`config/next-scan-v2.json`](config/next-scan-v2.json).
It uses micrOMEGAs 7.1.4, the observed LZ WS2024 numerical table and BSMPT 3.2.1.
All evaluated points are retained. DM, global-vacuum and RG results define
independent subsets; they do not gate EWPT exploration. Baryogenesis and GW
candidate flags follow the separate criteria documented there. Older workflow
examples below retain historical conventions where they differ from v2.

For micrOMEGAs-only work, start with the [DM installation instructions](DM/README.md)
and [standalone scan/card/plot walkthrough](DM/steer_example/README.md).
The [model change audit](DM/models/changes-v2.md) explains the loop corrections
and the roles of `x1` and the retained four-gluon auxiliary fields.

The independent [Upsilon flavour constraint](docs/flavour-constraints.md) is
enabled for new scans. Full viability and MG5 eligibility require a flavour
pass; EWPT eligibility keeps its existing definition. To add flavour results
to a saved `vx=0` scan without rerunning DM, HiggsTools or EWPT:

```bash
python reevaluate_trsm_flavour.py old_scan.dat --output scan_with_flavour.tsv
python plot_trsm_constraint_suite.py scan_with_flavour.tsv --format both
```

The input remains unchanged. Older files without flavour assessments are
labelled unassessed and cannot satisfy the updated full-viability selection.
Start a new scan after upgrading; pre-flavour campaign checkpoints cannot be
resumed with the changed physics/output prescription.

# Instructions:

## Automated installation

Use the [runtime bootstrap](docs/runtime-bootstrap.md) to download pinned
software releases, build both micrOMEGAs backends, HiggsTools and BSMPT, fetch
the Higgs datasets, and generate/compile all supported MadGraph processes:

```bash
python3.13 tools/bootstrap_runtime.py --prefix ../runtime-v2-new --dry-run
python3.13 tools/bootstrap_runtime.py --prefix ../runtime-v2-new --jobs 4
source ../runtime-v2-new/activate.sh
```

The installer requires a fresh prefix and supports `--resume`. Existing
runtimes and campaigns are preserved. See the linked guide for prerequisites,
component selection, provenance, and the MadGraph central-PDF prescription.
For an older MG5 installation with a missing `ml5_..._polynomial_constants.mod`
file, use the [MadLoop repair command](docs/runtime-bootstrap.md#repairing-an-existing-madloop-installation)
to rebuild the initialization artifacts while preserving saved campaigns.
For an older campaign failing with a zero MG5 survey or a `DiscreteSampler`
error, use the [survey repair instructions](docs/runtime-bootstrap.md#repairing-an-existing-mg5-survey).

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
fully viable points: the v2 tree-theory, assessed experimental, independent
`flavour`, and aggregate `dm` verdicts must all be `True`.

To drop only the DM requirement while retaining the theory, experimental, and
flavour constraints, use:

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

### Add loop-induced ZH2 rates to a saved campaign

The additional process `gg_eta0Z` evaluates `g g > eta0 z [noborn=QCD]`.
For new scans, request it explicitly alongside the existing channels using
`--mg5-process gg_heta0 --mg5-process pp_eta0Z --mg5-process gg_eta0Z`.
The default two-process scan selection remains unchanged.

For the saved Odysseus campaign, activate its usual runtime and run these
commands from the repository checkout containing the updated scripts:

```bash
source ../runtime-v2-new/activate.sh &&
python tools/setup_mg5_process.py --mg5-process gg_eta0Z &&

python reprocess_trsm_mg5.py \
  output/odysseus-100x1000-66662/combined_points.tsv \
  --output output/odysseus-100x1000-66662/combined_points_mg5.tsv \
  --energy 13.6 &&

python plot_trsm_constraint_suite.py \
  output/odysseus-100x1000-66662/combined_points_mg5.tsv \
  --output-dir plots/odysseus-100x1000-66662-updated \
  --format both
```

The setup command uses `TRSM_MG5_LOCATION` and adds only the requested process.
It validates an existing matching installation and refuses a conflicting one.
The reprocessor streams every stored row to a separate TSV, retains DM, BSMPT
and other assessments, and calculates missing rates only for the existing
version-aware full-viability selection. Finite existing production and signal
rates are reused. Add `--dry-run` to inspect the number of missing rates; after
an interruption, repeat the reprocessing command with `--resume` to continue
from its per-process SQLite checkpoint. The input, settings and runtime must
still match. This augments saved points; it does not resume the original scan
with a changed process list. If the source has no energy metadata, `--energy`
is the explicit assumption for the existing rates as well as the new rates.
The enriched TSV is published only after every requested rate finishes. The
`&&` operators stop the sequence on an earlier error, so plotting starts only
after its input exists. Use `--resume` only when a `.partial` checkpoint exists.

`mg5_xsec_pp_eta0Z_pb` and `mono_z_xsec_pb` retain their tree-level meanings.
The new columns are:

| Column | Meaning |
| --- | --- |
| `mg5_xsec_gg_eta0Z_pb` | Loop-induced gg production cross section |
| `mg5_xsec_pp_eta0Z_total_pb` | Tree production + loop-induced gg production |
| `mono_z_gg_xsec_pb` | gg production times stored `h2_h3h3_br` |
| `mono_z_total_xsec_pb` | Total production times stored `h2_h3h3_br` |

All rates are in pb. Both production components must be available for a total;
missing rates are not replaced by zero. The sum is labelled “tree +
loop-induced gg”, rather than being presented as a complete NLO prediction.

### Combined experiment + DM portal plot and HH enhancement

The additional dashboard
`dashboard_portal_mass_plane_summary_experimental_dm.png` shows `K133` and
`K233` for the intersection of the existing experimental and aggregate DM
selections. The original dashboard with separate selections remains available.
The new panels use the same coupling colour scales and gray background points
as the original portal maps.

`dashboard_hh_enhancement.png` shows the full-viability points on the mass
plane and versus `M2` and `M3`, with

```text
R_HH = mg5_xsec_hh_pb / 0.01617
```

This is the full LO `gg -> H1 H1` production rate before Higgs decays, divided
by the project's existing 13.6 TeV SM reference of 16.17 fb (MG5, `nn23lo1`).
The unit line denotes the SM rate. Use `--hh-sm-xsec-pb VALUE` to supply a
different matched SM reference at 13.6 TeV. An explicit incompatible energy
in the source metadata makes the HH plots unavailable; an input without
energy metadata is assumed to be at 13.6 TeV. Negative or nonfinite rates are
unavailable, and zero rates are counted but omitted from logarithmic axes.
The index and summary record the normalization, coverage, enhancement range
and number of points above the SM rate.

Both Odysseus files `combined_points.tsv` and `combined_points_mg5.tsv` were
checked on 21 October 2026: neither contains `mg5_xsec_hh_pb`. They contain
`xs136_lo_h2_pb` and `h2_h1h1_br`, which supply the factorized resonant
contribution. That contribution alone does not include the continuum or its
interference with the resonance and cannot give the full HH enhancement.
The enriched file has 591,340 stored rows, 835 passing experiment + DM, and
571 full-viability points requiring HH calculations.

From the updated repository checkout on Odysseus, run:

```bash
source ../runtime-v2-new/activate.sh &&
python tools/setup_mg5_process.py --mg5-process hh &&
python reprocess_trsm_mg5.py \
  output/odysseus-100x1000-66662/combined_points_mg5.tsv \
  --output output/odysseus-100x1000-66662/combined_points_mg5_hh.tsv \
  --mg5-process hh --energy 13.6 &&
python plot_trsm_constraint_suite.py \
  output/odysseus-100x1000-66662/combined_points_mg5_hh.tsv \
  --output-dir plots/odysseus-100x1000-66662-portal-hh \
  --format both \
  --plots dashboard_portal_mass_plane_summary_experimental_dm dashboard_hh_enhancement
```

The reprocessor calculates only missing HH rates for full-viability points
and writes a separate TSV. Existing rates and constraint assessments are
preserved. After an interruption, repeat its command with `--resume`; completed
HH calculations are read from the checkpoint. To inspect the workload first,
add `--dry-run` to the reprocessing command. Its default process list remains
the three associated-production channels; HH must be requested explicitly.

The `--plots` option generates only the named dashboards and their individual
panels, together with an index and summary. Omit it to generate the full
suite. To make the combined portal plot immediately, use the existing
`combined_points_mg5.tsv` with just
`--plots dashboard_portal_mass_plane_summary_experimental_dm`.
PDF versions of both dashboards are also produced with `--format both`.

### Determinant versus transition strength and CalcGW

`dashboard_bsmpt_determinant` adds two plots with dimensionless `D` on the
x axis: the EW-entry jump `Delta v_EW(Tc)/Tc`, and the maximum any-field jump
`Delta phi/T` over the stored critical, nucleation and percolation diagnostics.
The first is the requested EWPT strength plot. Completion is not required for
this critical-temperature diagnostic. The `_full_viability` dashboard applies
the existing theory, experiment, flavour and aggregate DM selection as well.
Both signs of D and zero are displayed using a signed logarithmic axis with
linear threshold `1e-5`; no determinant cut is added.

D is calculated from the saved parameters using the
`test_trsm_ewpt.check_eq_4_18` convention, matching the quartic matrix in
[the analytic EWPT treatment](https://arxiv.org/html/2404.12425v2):

```text
D = c_phi c_X c_S + 2 c_phiX c_phiS c_XS
    - c_phiS^2 c_X - c_phi c_XS^2 - c_phiX^2 c_S
c_ii = lambda_i
c_ij = lambda_ij / 2   (i != j)
lambda_phi = (M1^2 cos(a12)^2 + M2^2 sin(a12)^2) / (2 v^2)
lambda_S   = (M2^2 cos(a12)^2 + M1^2 sin(a12)^2) / (2 vs^2)
lambda_phiS = (M2^2 - M1^2) cos(a12) sin(a12) / (v vs)
```

The other three quartics are stored as `lX`, `lPhiX`, `lSX`; `M1` and `v`
come from the shared SM inputs. This is the named Eq. (4.18) diagnostic, not
the determinant of the scalar mass matrix.

The normalization was corrected on 29 October 2026: earlier code also halved
the diagonal entries. The correction applies to both plotters and the optional
`--ewpt-require-eq418` check, including its principal minors. It can change the
sign of D and the prefilter decision. Regenerate D plots from the saved scan
parameters; existing images and saved prefilter outcomes are not rewritten.
Points previously skipped by that optional filter need reassessment if its
updated selection is required. With both X portals zero, the corrected
determinant satisfies `D = lX * M1^2 * M2^2 / (4 v^2 vs^2)` for any `a12`.

The portal dashboard was checked against all 835 experiment + DM points in
the Odysseus campaign on 22 October 2026: its left panel uses `K133`, and its
right panel uses `K233`, with no interchange. The columns exactly reproduce
`K133 = (lPhiX*v*cos(a12) - lSX*vs*sin(a12))/2` and
`K233 = (lPhiX*v*sin(a12) + lSX*vs*cos(a12))/2`. These are potential
coefficients in GeV for `h1*h3^2` and `h2*h3^2`, respectively; the vertices
are `-2i K133` and `-2i K233`. Both dashboard axes are masses, not couplings.

The HH-enriched campaign now exists as `combined_points_mg5_hh.tsv`. It has
591,340 stored rows and 31,204 points with `ewpt_eligible=True`, but no native
CalcGW output columns. `reprocess_trsm_gw.py` runs BSMPT **3.2.1 CalcGW** on
those saved eligible points, including points without a prior strong-transition
flag. It uses the same upper tracing temperature, 1000 GeV. CalcGW recomputes
the thermal and bounce calculation from the parameters; it cannot consume the
saved CalcTemps results as a bounce cache. No new random scan is generated.

The defaults are the native percolation choice `T*=Tp`, wall speed `vwall=0.95`
and turbulence parameter `epsturb=0.1`. These can be changed with `--trans-temp`,
`--vwall` and `--epsturb` for a separate output. Native percolation/completion
probabilities, minimizers and other numerical settings keep the installed
BSMPT defaults. Alpha and beta/H are the native quantities:

```text
alpha = [Delta V - (T*/4) Delta(dV/dT)] / rho_radiation(T*)
beta/H* = T* d(S3/T)/dT evaluated at T*
Delta = false vacuum minus true vacuum
```

Each transition has its own `gw_transitions` JSON entry, including native
status, temperature, alpha, beta/H, efficiency factors, characteristic
frequencies/amplitudes and LISA SNR (three years). Raw TSVs, logs and commands
are saved in `<output>.work/row-*`. `gw_native_summary` also records the
thermal statuses and transition history. A plotted transition must have
successful native GW and requested-temperature statuses, successful NLO,
tracing and coexistence assessments, finite positive alpha and T*, and
beta/H >= 1. Failed/time-out calculations are recorded as such, not zero
signals or physical exclusions. All old scan cells and metadata are copied
into a separate output, with the new configuration and executable/source
hashes recorded in its metadata.

`dashboard_calcgw` shows alpha versus beta/H, each versus M2, the sound-wave
second characteristic frequency versus amplitude, LISA SNR versus M2, and
alpha versus T*. Each marker represents one successful coexisting-phase
transition; multiple markers can belong to a point. No selection of a
cosmological transition history is imposed. The sound-wave panel is a
component diagnostic, not a total GW spectrum. The `_full_viability` version
applies the stored full-viability selection. Nonpositive/nonfinite values
are omitted from logarithmic panels and counted. The original `ewpt_gw_candidate`
strength screen remains separate from these native predictions. See the
[BSMPT v3 paper](https://arxiv.org/abs/2404.19037) for the thermal/GW framework;
the adapter uses the installed 3.2.1 output schema and defaults.

All CalcGW panels now use stars for the same parent scan points as plot 34b:
successful stored critical EW entry with `ewpt_ew_entry_jump_over_T > 1`.
Every plotted transition from such a point gets a star; that particular GW
transition need not itself be EW entry. Missing diagnostics do not create stars.
The marker colours continue to show temperature (or beta/H), including in
the full-viability panels and dashboards.

Plot 97 includes a horizontal **SNR = 5** reference for the native three-year
LISA prediction. This is a sensitivity threshold, not a new viability cut or
a five-sigma discovery claim. New plot
`97b_calcgw_lisa_snr_vs_m2_1_to_100` displays the same data with the SNR axis
limited to 1--100; it also has a `_full_viability` version. Values outside
the window are not moved to its edges. The annotation counts the transitions
and distinct points in view.

For the saved `odysseus-100x1000-66662` GW enrichment, **plot 97 without full
viability already includes the experimental constraints**: all 20,763 points
with successful CalcGW results explicitly pass tree theory, HiggsBounds,
HiggsSignals, EWPO, the applicable W-mass comparison and flavour. The general
EWPT gate permits unknown experimental assessments, but this campaign has none
among its eligible points. The full-viability selection leaves 431 successful
GW points (647 transitions), by additionally requiring an aggregate DM pass.
No saved transition reaches SNR = 5. The largest SNR is 1.5344, at
`M2 = 401.3081 GeV`, `M3 = 62.5575 GeV`, on a starred point that fails DM.
It is the only transition in the 1--100 window; the full-viability window is
empty, with a maximum SNR of `4.4912e-11`.

Regenerate the six CalcGW panels, both dashboards and the two zooms from the
complete saved ledger with:

```bash
python plot_trsm_constraint_suite.py \
  output/odysseus-100x1000-66662/combined_points_mg5_hh_gw.tsv \
  --output-dir plots/odysseus-100x1000-66662-calcgw-20261029 --format both \
  --plots dashboard_calcgw dashboard_calcgw_full_viability \
  97b_calcgw_lisa_snr_vs_m2_1_to_100 \
  97b_calcgw_lisa_snr_vs_m2_1_to_100_full_viability
```

The collaborator constraint checklist is [applied-constraints.tex](docs/applied-constraints.tex).
Its included [HiggsTools tables](docs/higgstools-analyses.tex) list the 31
HiggsBounds IDs with nonzero selection counts and all 22 HiggsSignals measurement files, with
scope, process, coverage and references. The selected-point counts use all
591,340 reference-campaign rows; 31 different HB IDs were selected somewhere
in that sample. Rows with zero selection counts are omitted from the tables;
the JSON inventory includes all 310 installed HiggsBounds IDs. The tables do
not treat all observed bounds as simultaneous cuts. To rebuild the tables and PDF:

```bash
python tools/export_higgstools_tables.py \
  --inventory validation/collaborator-update-20261029/higgstools-inventory.json \
  --selected-counts validation/collaborator-update-20261029/higgs-selected-counts.json \
  --selected-only
cd docs
pdflatex applied-constraints.tex
pdflatex applied-constraints.tex
```

Omit `--selected-only` to include zero-count rows; omit both `--selected-only`
and `--selected-counts` to export the installed catalogue without historical
selection counts. The two new mass-range profiles and their 100-worker launch
commands are documented in [the parallel campaign guide](docs/parallel-campaigns.md#mass-range-runs-from-the-66662-reference).

The checked [PDF](output/pdf/applied-constraints.pdf) and
[standalone LaTeX](output/pdf/applied-constraints.tex) are available directly.
The [collaborator handover kit](output/odysseus-66662-handover-kit.zip) includes
these documents, detailed scan metadata, a Python reader and
[reading instructions](output/collaborator-66662/READ_ME.txt). Send it with the
final `combined_points_mg5_hh_gw.tsv`, its metadata sidecar and
`campaign_state.json`. The reader uses only the Python standard library and
exports a small overview CSV plus complete rows for the fully viable points:

```bash
python output/collaborator-66662/inspect_points.py \
  output/odysseus-100x1000-66662/combined_points_mg5_hh_gw.tsv \
  --out output/quicklook-66662
```

Raw campaign files, local plot directories and native scratch output are
ignored by Git. The repository includes the compact handover documents and
validation records needed to interpret the published changes.

After updating this checkout on Odysseus, make the determinant and portal
figures immediately:

```bash
cd /home/apapaefs/Projects/TwoSingletDM/twosingletDM
source ../runtime-v2-new/activate.sh &&
python plot_trsm_constraint_suite.py \
  output/odysseus-100x1000-66662/combined_points_mg5_hh.tsv \
  --output-dir plots/odysseus-100x1000-66662-determinant \
  --format both --plots dashboard_bsmpt_determinant \
  dashboard_bsmpt_determinant_full_viability \
  dashboard_portal_mass_plane_summary_experimental_dm
```

Build the missing CalcGW target in the existing runtime, run a three-point
pilot, then resume across the eligible campaign points:

```bash
cmake --build ../runtime-v2-new/BSMPT-3.2.1/build --target CalcGW -j 8 &&
python reprocess_trsm_gw.py \
  output/odysseus-100x1000-66662/combined_points_mg5_hh.tsv \
  --output output/odysseus-100x1000-66662/combined_points_mg5_hh_gw.tsv \
  --thigh 1000 --jobs 1 --max-points 3 &&
python reprocess_trsm_gw.py \
  output/odysseus-100x1000-66662/combined_points_mg5_hh.tsv \
  --output output/odysseus-100x1000-66662/combined_points_mg5_hh_gw.tsv \
  --thigh 1000 --jobs 8 --resume &&
python plot_trsm_constraint_suite.py \
  output/odysseus-100x1000-66662/combined_points_mg5_hh_gw.tsv \
  --output-dir plots/odysseus-100x1000-66662-with-gw --format both
```

Each CalcGW process is single threaded; `--jobs` bounds concurrent processes.
`--timeout` defaults to 3600 seconds per point. `--dry-run` only counts the
workload. `--max-points` limits new calculations for a pilot; an incomplete
pilot writes a checkpoint, without publishing a partial final TSV. Repeat
with `--resume` and without the limit to finish. Each result, including
native failures and timeouts, is checkpointed; resume does not retry completed
attempts. A changed timeout or physics configuration requires a new output.
Existing input, output and checkpoint identities are checked before reuse.
The final command generates the entire suite, including HH and the new
figures. To plot only the additions, append `--plots dashboard_bsmpt_determinant
dashboard_bsmpt_determinant_full_viability dashboard_calcgw
dashboard_calcgw_full_viability dashboard_portal_mass_plane_summary_experimental_dm`.

The gravitational-wave FOPT candidate-status legend sits below the plotting
area, including in the BSMPT dashboards and full-viability versions. To
regenerate these figures from the saved GW-enriched results, run:

```bash
python plot_trsm_constraint_suite.py \
  output/odysseus-100x1000-66662/combined_points_mg5_hh_gw.tsv \
  --output-dir plots/odysseus-100x1000-66662-bsmpt-layout-fixed \
  --format both --plots dashboard_bsmpt_summary dashboard_bsmpt_summary_full_viability
```

This also regenerates the individual panels. The separate output directory
keeps the existing full-suite index available. No CalcGW rerun is needed.

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

Select the micrOMEGAs backend with `--micromegas-version 7` (7.1.4) or
`--micromegas-version 6` (6.1.15). The exact version strings also work.
New v2 commands default to 7.1.4. For example, from this repository:

```bash
python generate_trsm_points.py 123 --nrandom 500 --micromegas-version 7
```

Each backend uses `../runtime-v2/micromegas_<version>/TRSM/main` relative to the
repository. Set `TRSM_RUNTIME_ROOT` to select a different runtime parent.
For another installation, add `--micromegas-main /absolute/path/to/TRSM/main`
and specify its version with `--micromegas-version`. The executable is checked
before any scan output is created. The version and resolved executable path
are recorded in the metadata sidecar. Version 7 run tags include `-mo7.1.4`,
so comparisons using the same date and seed have separate output files.
Custom executable paths add `-customMO`.

New version 7 scans enable the Planck CMB annihilation constraint by default
and add `-cmb-planck2018` to the filename. Use `--no-planck-cmb` to disable
it. Version 6 is an explicit comparison backend, with CMB disabled; explicit
`--planck-cmb` requires a rebuilt, capable driver. The constraint uses the
built-in low-velocity s-wave approximation and rescales its bound ratio by
`min(1, Omega_h2 / 0.12)^2`. The existing relic-density upper cut remains
0.121. An unavailable CMB result prevents aggregate DM acceptance while
retaining valid relic-density and detection diagnostics. See
[the CMB treatment and validation](DM/planck-cmb.md).

Resume with `--resume-from` alone: the saved backend is restored and cannot
be changed during a campaign. Older metadata/checkpoints without a backend
selection retain their implicit 6.1.15 backend and remain compatible.
Campaigns without saved CMB settings resume with CMB disabled, including
older version 7 campaigns. Their column layouts and fingerprints are preserved.
CMB settings cannot change during resume.

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
fraction arrays remain the pre-invisible/base inputs; the HiggsTools adapter
combines them with the true exotic partial widths before setting the physical
branching fractions.

HiggsTools internally treats positive total widths below `1e-10 GeV` as zero
and otherwise rejects branching-ratio input for those points. Such widths are
therefore floored to `1e-10 GeV` only in the HiggsTools particle record. The
physical `w1`/`w2` values written to scan output and MadGraph parameter cards
retain their unfloored values.

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

When `--run-ewpt` is enabled, the viable-point TSV includes the original
selected-transition diagnostics and separate qualitative candidate flags:

```text
ewpt_ew_true_over_T
ewpt_ew_jump_over_T
ewpt_ew_entry_true_over_T
ewpt_ew_entry_false_over_T
ewpt_ew_entry_jump_over_T
ewpt_ew_entry_temperature_GeV
ewpt_ew_entry_temperature_kind
ewpt_ew_entry_transition_index
ewpt_ew_entry_nucl_jump_over_T
ewpt_ew_entry_nucl_temperature_GeV
ewpt_ew_entry_perc_jump_over_T
ewpt_ew_entry_perc_temperature_GeV
ewpt_ew_entry_percolated
ewpt_ew_entry_completed
ewpt_baryo_candidate
ewpt_gw_crit_field_jump_over_T
ewpt_gw_crit_temperature_GeV
ewpt_gw_nucl_field_jump_over_T
ewpt_gw_nucl_temperature_GeV
ewpt_gw_perc_field_jump_over_T
ewpt_gw_perc_temperature_GeV
ewpt_gw_max_field_jump_over_T
ewpt_gw_max_temperature_kind
ewpt_gw_max_transition_index
ewpt_gw_candidate
```

The first two values describe the same selected BSMPT transition. It is chosen
from the available strengths with priority `nucl`, then `perc`, then `compl`, then
`crit`; within one temperature kind the largest finite `ew_true/T` is selected.
If EWPT is not run, Eq. 4.18 skips the run, or no finite strength is available,
the columns are written as `nan`.

For the baryogenesis candidate flag, a critical-temperature CalcTemps FOPT
must cross from `|w1_false| < --ewpt-w1-threshold` to
`|w1_true| >= --ewpt-w1-threshold` (default 5 GeV). If several transitions
do so, the scan reports the one with the largest
`|w1_true-w1_false|/T_c`. `ewpt_baryo_candidate` is `True` precisely when
this critical-temperature EW jump is **strictly greater than 1**. It does
not require a nucleation, percolation, or completion result. The same
transition's EW jump and temperature at nucleation and percolation are
reported separately when BSMPT supplies them; they never change the
baryogenesis flag. The true- and false-phase EW VEV ratios and the
percolation/completion indicators are diagnostics only.

For the separate gravitational-wave candidate flag, define the total
field-space VEV jump
`Delta phi = sqrt((Delta w1)^2 + (Delta wx)^2 + (Delta ws)^2)`.
For **any** CalcTemps FOPT, including a singlet-only transition starting
from zero VEV, `ewpt_gw_candidate` is `True` if
`Delta phi/T > 1` at the critical, nucleation, **or** percolation
temperature. The three temperature-specific columns give the largest
ratio over all FOPTs at each temperature, so they can refer to different
transitions. The maximum ratio, its temperature kind, and its CalcTemps
transition index are also stored. Approximate nucleation and completion
results do not enter either candidate cut.

Both flags are qualitative scan labels, not baryogenesis or gravitational-wave
predictions. In particular, a field jump alone does not establish a
gravitational-wave signal, and a large EW jump alone does not calculate the
sphaleron washout rate. BSMPT nucleation/percolation estimates are retained
for inspection without being required for the baryogenesis candidate flag.
These new flags are additional annotations; the existing viability columns
and the selected-transition diagnostics keep their previous definitions.

The scan also saves micrOMEGAs `Xf = m_DM/T_f` as `dm_xf` and the derived
`dm_freezeout_temperature_GeV = dm_mdm/dm_xf`. A failed or unavailable DM
calculation leaves both as `nan`; no approximate `m_DM/20` value is inserted.
These are the freeze-out diagnostics from the standard unbroken-\(Z_2\)
micrOMEGAs calculation. The [micrOMEGAs manual](https://lapth.cnrs.fr/micromegas/v6.0/manual_6.0.pdf)
defines `Xf` through the abundance departure from equilibrium.
To populate these fields in an older scan, rerun the DM calculation with
`reevaluate_trsm_dm_higgs.py`; an EWPT-only reprocessing run cannot infer a
missing `Xf` from the existing TSV.
Re-evaluating DM clears any previously stored freeze-out/phase comparison
flags and thermal-VEV diagnostics, since a changed `Xf` would make them stale; a subsequent EWPT run
with `--rerun-existing-ewpt` recomputes them.

For points with a MinimaTracer global branch, the scan records the lowest and
highest *sampled* temperatures of an \(X\)-broken global minimum in
`ewpt_x_broken_min_T_GeV` and `ewpt_x_broken_max_T_GeV`. The JSON-valued
`ewpt_x_broken_intervals_GeV` retains separate sampled windows if the history
has more than one. `ewpt_x_final_restoration_low_T_GeV` and
`ewpt_x_final_restoration_high_T_GeV` bracket the final sampled return to
`wx=0` on cooling when the trace reaches an unbroken `T=0` vacuum. They are
`nan` if no such return is resolved. `ewpt_x_phase_at_freezeout` is `broken`, `unbroken`,
`boundary_unresolved`, or `outside_traced_range` when `T_f` is available.
`ewpt_x_broken_at_or_after_freezeout=True` means at least one global
\(X\)-broken sample lies at or below `T_f`; `False` requires a resolved
unbroken phase at `T_f`, coverage to \(T=0\), and no such sample. An
unresolved comparison is `nan`. The inverse flag,
`dm_relic_z2_freezeout_compatible`, is `True` when the sampled equilibrium
history has restored \(Z_2\) by nominal freeze-out, `False` when an \(X\)-broken
global sample lies at or below freeze-out, and `nan` when the comparison is
unresolved. It is a compatibility diagnostic, not proof that the micrOMEGAs
relic abundance is valid: the actual transition may differ from the
equilibrium global-minimum path, and re-equilibration after an earlier
\(X\)-broken phase is not checked. Local \(X\)-broken minima are excluded from
this thermal-history diagnostic. The plot
`35g_freezeout_vs_x_broken_window` compares `T_f` with the sampled lower and
upper window temperatures. These diagnostics do not change the DM, EWPT,
baryogenesis, or GW selections. A temporary \(X\)-broken phase can make the
standard relic calculation unreliable; the sampled comparison alone does not
recalculate the abundance through that phase.

The separate fixed-vacuum input screen compares the magnitudes of the global
branch's electroweak (`w1`) and singlet (`ws`) VEVs with their sampled (T=0)
values. It records their ratios at nominal `T_f` and the largest fractional
shifts in the window `[T_f/2, 2*T_f]`. The Boolean
`dm_relic_thermal_vev_shift_ge_10pct` is `True` when either observed shift is
at least 10%, `False` only with complete same-phase window coverage and both
shifts below 10%, and `nan` otherwise. The separate
`dm_relic_thermal_phase_boundary_bracket_overlaps_window` identifies a change
of global minimum whose sampled temperature bracket overlaps that window;
interpolation never crosses such a bracket. This 10% threshold and the window
are screening conventions, not relic-density validity conditions.

For every `vx=0` point the scan also records the signed zero-temperature gaps
`dm_resonance_h1_mass_gap_GeV = M1 - 2*M3` and
`dm_resonance_h2_mass_gap_GeV = M2 - 2*M3`, their absolute values divided by
`T_f` when it is available, and the mediator with the smaller absolute gap.
Positive gaps put the pole above the two-DM threshold. These ratios do not
include mediator widths or thermal mass shifts and are not resonance
significance tests. A small VEV shift can still matter near a narrow pole or
when the terms in `K133` nearly cancel. Neither this screen nor the $Z_2$
diagnostic changes the stored DM selection.

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

An explicitly selected `--dm-limit-table /path/to/verified-si-limits.json`
replaces the inherited `lz2025-source` direct-detection fit with a published
90% observed elastic SI upper-limit table. This works with both micrOMEGAs
versions and preserves the abundance rescaling. Tables are checked for units,
mass coverage and provenance; their checksum is recorded for checkpoint/resume.
See [the table format and LZ 2026 status](DM/direct-detection.md).
The official mass limits from arXiv:2609.02823 remain unavailable: the linked
HEPData record denied access even after login on 2026-09-17. An optional
[temporary 400-4000 GeV approximation](DM/data/lz2026/README.md) uses the
digitized Figure S7 elastic limit at 1 TeV and assumes a limit proportional
to DM mass. To select it, add
`--dm-limit-table DM/data/lz2026/lz2026-figs7-highmass-approx.json`
and restrict a random scan with `--m3-min 400 --m3-max 4000`.
It is explicitly labelled `highmass-approx`; the existing default is unchanged.

`--run-ewpt-on-dm-failed` is an exploratory option for otherwise-good points
that fail only the dark-matter check. It runs BSMPT for those points and writes
them to the `_dm_failed` sidecar, including `ewpt_ew_true_over_T` and
`ewpt_ew_jump_over_T` and the two candidate flags when BSMPT returns
transition data. Use it together
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
TSV. The EWPT diagnostics and candidate columns are added or updated only for
selected rows.
Rows already containing an EWPT attempt are preserved by default, so a scan
originally run with `--run-ewpt` does not repeat its fully viable BSMPT jobs;
add `--rerun-existing-ewpt` when those rows need the new candidate columns.
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

A separate gradual-DM dashboard starts with the existing **non-DM viability**
selection (version-aware theory, experimental and flavour passes), then applies
relic density, direct detection, indirect detection and CMB in that order.
It shows the `M2`–`M3`, `M3`–`K133` and `M3`–`K233` planes. Missing component
assessments cannot pass; an assessed indirect-detection result outside coverage
retains its non-veto meaning, and disabled CMB leaves the preceding population
unchanged. Counts and discrepancies with stored full viability are reported in
the summary. The original collaborator cumulative sequence remains separate.

When the scan contains recorded BSMPT results, plots 30--36 and a fifth
dashboard are added automatically. New scans add separate baryogenesis and
gravitational-wave candidate maps and temperature diagnostics:

- BSMPT run/failed/no-selected-FOPT/weak-FOPT/strong-FOPT status on the
  \(M_2,M_3\) plane;
- EW-entry \(\Delta v_{\rm EW}(T_c)/T_c\) and the baryogenesis candidate
  flag on the \(M_2,M_3\) plane (31c, 31d);
- the largest any-field \(\Delta\phi/T\) and the gravitational-wave FOPT
  candidate flag on the same mass plane (31e, 31f);
- the original selected \(v_{\rm EW,true}(T_*)/T_*\) and
  \(\Delta v_{\rm EW}(T_*)/T_*\) diagnostics (31, 31b);
- MinimaTracer global phase-route and electroweak-entry-step maps;
- the EW-entry critical jump versus \(M_2\) and \(M_3\) (34b, 35b),
  and the largest any-field jump versus \(M_3\) (35d);
- selected true EW VEV ratio versus the EW-entry critical jump (35c);
- nucleation and percolation jumps versus critical jumps, for the EW-entry
  transition and for any-field FOPTs separately (35e, 35f); and
- BSMPT evaluation, phase-history, and transition-result counts.

The generator selects \(T_*\) with priority nucleation, then percolation,
completion, and critical temperature. That selected transition may be a later
broken-to-broken step, so its \(v_{\rm EW,true}(T_*)/T_*\geq1\) split must not
be read as evidence that EW symmetry was strongly broken when it first appeared.
The independent baryogenesis and gravitational-wave candidate counts are
recorded in `constraint_summary.tsv` when their columns are present. The
nucleation and percolation estimates are plotted for comparison and do not
gate the baryogenesis flag.
The phase-route map distinguishes direct electroweak entry, paths through an
\(S\)-broken phase, paths through an \(X\)-broken phase, and other multistep
histories. Since an \(X\)-broken phase violates the nominal dark-sector
\(Z_2\), those points are also counted explicitly in `constraint_summary.tsv`.
Rows for which BSMPT was not requested remain visible as “not run”, while
failed evaluations are never interpreted as successful or as having no
first-order transition.

For a file with no recorded BSMPT attempt, these figures and their
dashboard are skipped. Cascade/MadGraph figures are independently skipped when
their named scan columns are absent, so legacy inputs remain usable; the HTML
index explains every unavailable plot. A finite `ewpt_ew_true_over_T` is
sufficient for compatibility with older files. Newer files additionally use
`ewpt_ew_jump_over_T`, `ewpt_status`, `ewpt_global_phase_path`,
`ewpt_has_x_broken`, and `ewpt_ew_step_index`. Candidate plots are omitted
for legacy scans without their respective columns.
Detailed transition
temperatures remain in the per-point `ewpt_result.json` files and are not
reconstructed by the scan-table plot suite.

Across BSMPT scatter plots, **stars identify strong critical EW entry**:
successful assessment, `|vEW,false| < 5 GeV <= |vEW,true|`, and strictly
`Delta vEW/Tc > 1`, using the stored EW-entry diagnostic. Neither completion
nor percolation is an additional requirement. A large selected broken-to-broken
transition, a GW-only candidate, or a multistep history does not by itself earn
a star. Their other classifications and stored flags remain available. Legacy
inputs without EW-entry diagnostics explicitly show the star classification
as unavailable.

All 19 BSMPT plots and their dashboard have additional `_full_viability`
counterparts in the **BSMPT + DM full viability** section. Their points, counts,
legends and denominators use full viability. Empty or unavailable selections
are explained in the report instead of being filled with unselected points.

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

The full-viability mono-Higgs/Z plots appear under **Full-viability collider
signal plots**, independently of YR4 signal availability. Additional figures
show gg-only and total mono-Z rates versus both masses, and compare tree, gg
and total raw ZH2 production. Aggregate-DM variants of the K133/K233-versus-M3
resonance plots complement the existing experimental and relic-only variants.

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

For a full reevaluation with version 7 and its default CMB constraint:

```bash
trsmdm/bin/python reevaluate_trsm_dm_higgs.py \
  output/trsm_points_OLD.dat \
  --output output/trsm_points_OLD_mo7_cmb.dat \
  --micromegas-version 7 --checkpoint-every 25
```

The same `--planck-cmb` / `--no-planck-cmb` controls apply. Reevaluation
recovers the input's direct-detection treatment from its row provenance and
metadata. A tabulated limit must have a recoverable, matching source file;
otherwise supply `--dm-limit-table /path/to/limits.json` explicitly. This
also preserves the temporary LZ 2026 approximation when it was selected.
An intentional change of SI treatment can use `--dm-limit-table` or
`--dm-limit-model lz2025-source` / `--dm-limit-model legacy-output`.

Repeat the evaluation options when resuming. The checkpoint verifies the
backend, SI-table hash and CMB treatment; the completed file has a metadata
sidecar recording them. Legacy reevaluation checkpoints resume with CMB
disabled and the original default SI fit, then record the selected backend
for subsequent resumes. CMB fields in an input are always recomputed or
cleared and marked disabled. Older files are shown as CMB-unassessed by the
constraint plotting suite. No existing scan is reevaluated automatically.

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

The helper resolves the BSMPT build relative to the scan repository, independent
of the working directory:

```text
../BSMPT/build/macos-armv8-release/bin/CalcTemps
../BSMPT/build/macos-armv8-release/bin/MinimaTracer
```

On manto this is `/Users/apapaefs/Projects/BSMPT/build/macos-armv8-release/bin`.
On the laptop it is
`/Users/apapaefs/Projects/TwoSingletDM/BSMPT/build/macos-armv8-release/bin`.
These defaults also apply to `generate_trsm_points.py` and
`reprocess_trsm_ewpt.py`, including the EWPT step enabled by
`--write-all-points`. No executable flags are needed for either layout.

Override the binaries if needed (use `--ewpt-executable` and
`--ewpt-minima-executable` in the scan and reprocessing scripts):

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

and inspect the critical EW-entry jump `ew_jump/T`. The scan's qualitative
baryogenesis flag uses this ratio at `T_crit` and does not require a
successful bounce, nucleation, or completion calculation. The latter
results are still saved for physical interpretation and possible
gravitational-wave follow-up.

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

For configurable concurrency, evo/thc targets, automatic campaign resume and
Odysseus launch/pilot commands, use the [parallel campaign guide](docs/parallel-campaigns.md).
Both campaign launchers support `--run-mg5`, repeated `--mg5-process` selections,
and `--mg5-without-dm`. MG5 runs use one core and an exclusive lock per shared
process directory. Use a fresh campaign to change its MG5 settings.
The examples below retain the older basic CLI; `--run-cwd` now selects a parent
for isolated per-seed directories.

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


## Install micrOMEGAs 6.1.15 or 7.1.4 and create TRSM

Both releases use the same tracked `DM/main.c` driver, canonical `h4GOn`
model, card mapping, standard-cosmology `darkOmega` calculation, and Python
relic/direct/indirect exclusion treatment. The driver keeps `VWdecay=0` and
`VZdecay=0`, signed nucleon amplitudes, and the Fermi-LAT R16 photon-line
calculation. Selecting version 7 does not activate its optional nonstandard
cosmology or additional experimental constraints. Predictions can change with
the upstream release even with these settings held fixed.

Download 7.1.4 from the [official installation page](https://micromegasdm.github.io/v7.1/install_7.1.html)
and extract it beside this repository. From the repository root:

```bash
curl -fL https://micromegasdm.github.io/downloadarea/v7.1/micromegas_7.1.4.tgz \
  -o ../micromegas_7.1.4.tgz
tar -xzf ../micromegas_7.1.4.tgz -C ..
sh DM/setup_micromegas.sh ../micromegas_7.1.4
../micromegas_7.1.4/TRSM/main ../micromegas_7.1.4/TRSM/data.par
```

The source archive used for the manto installation has SHA-256
`c8cf207b17541a5b7d7e7ff157f1c1eb36e1dd3f7547e7745559b195c2a34264`.
The same setup script supports a fresh 6.1.15 extraction. It builds serially
(upstream Makefiles write shared static archives), creates `TRSM`, installs
the model and driver, and applies the version-specific finite-value guard in
`DM/patches`. The guard returns an error instead of looping in channel sorting
when failed calculations leave non-finite annihilation weights. It is still
needed in 7.1.4. The script refuses to overwrite an existing `TRSM` directory.

On manto the executables are
`/Users/apapaefs/Projects/micromegas_6.1.15/TRSM/main` and
`/Users/apapaefs/Projects/micromegas_7.1.4/TRSM/main`.
See [the manto validation record](DM/micromegas7-manto-validation.md) for the
benchmark comparison, test results, and retained log locations.

The checked-in `h4GOn` and `h4GOff` models implement the canonical-v1
normalization: the potential contains `LHX/2`, `LSX/2`, and `LX/4`, while the
input card remains one-to-one (`LHX=lPhiX`, `LSX=lSX`, and `LX=lX`). Do not
compensate by scaling card values. Back up an existing model and executable
before replacing them, and do not rebuild an installation in use by a scan.

## run micrOMEGAs with a steering code:
- check out DM/example_steer/README.md and create another directory for the study you'd like to perform.
