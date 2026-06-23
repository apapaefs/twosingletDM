# twosingletDM
TRSM + Dark Matter + ElectroWeak Baryogenesis vs. Higgs Boson Pair production

# Instructions:

## Download MG5, e.g. 2.9.22 (https://launchpad.net/mg5amcnlo)
- Copy ```loop_sm_twoscalar_generic.tar.gz``` into MG5_aMC_2_9_22/models and untar it: ```tar xvzf loop_sm_twoscalar_generic.tar.gz```
- Launch MG5 and generate the process: 
```
./bin/mg5_aMC
import model loop_sm_twoscalar_generic
generate g g > h h [noborn=QCD]
output gg_hh_twoscalar
launch
```
- Enter and proceed to next screen, edit run card and change to the desired beam energies. -
- Exit and edit ```MG5_aMC_2_9_22/input/mg5_configuration.txt``, changing:
```
automatic_html_opening = False
```
- In the ```generate_mg5_trsm_xsecs.py``` script, change the following to the absolute directory of MG5:
```
MGLocation = '/home/apapaefs/Projects/TwoSingletDM/twosingletDM/MG5_aMC_v2_9_22/'
```
- Make sure ```ProcLocation['hh'] = 'gg_hh_twoscalar/'``` is set to the directory in which you have outputted the process.
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
K133_min = 1E-4
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
K133_pow_min = -3
K133_pow_max = 3

K233_pow_min = -3
K233_pow_max = 3
```

In either K-scan mode, each sampled `K133` and `K233` point is converted to the
potential-basis inputs expected by the existing BSMPT and micrOMEGAs pipeline:

```text
lphix = 2 * (cos(a12) * K133 + sin(a12) * K233) / 246
lsx   = 2 * (-sin(a12) * K133 + cos(a12) * K233) / vs
```

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
--resonantDM1: m3 = 2*m1 = 2*125.09 GeV
--resonantDM2: m3 = 2*m2
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

In a random scan, `--resonantDM1` fixes `m3` to the Higgs-resonant value for all
points, while `--resonantDM2` updates `m3` point-by-point after each random
`m2` draw.

To scan uniformly near either resonant mass relation instead of exactly on it,
use `--approximate-resonantDM` with a mass-window half-width:

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

The approximate mode is mutually exclusive with `--resonantDM1` and
`--resonantDM2`. The `--delta-res` value is in GeV.

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
```

This is selected from the available BSMPT strengths with priority `nucl`, then
`perc`, then `compl`, then `crit`. If EWPT is not run, Eq. 4.18 skips the run, or
no finite strength is available, the column is written as `nan`.

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
fail the dark-matter check.

`--run-ewpt-on-dm-failed` is an exploratory option for otherwise-good points
that fail only the dark-matter check. It runs BSMPT for those points and writes
them to the `_dm_failed` sidecar, including `ewpt_ew_true_over_T` when BSMPT
returns a finite strength. Use it together with `--run-ewpt` if you want BSMPT
for both viable and DM-failed points; by itself it only targets DM-failed points.

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
python3 -m py_compile generate_trsm_points.py run_trsm_seed_campaign.py
```


## Install micrOMEGAs 6.1.15 and create TRSM:

- ```cd DM``` and untar `micromegas_6.1.15.tar`: ```tar xvzf micromegas_6.1.15.tar```

- ```cd micromegas_6.1.15/``` then follow the READ.ME instructions I and V:
```[g]make``` and ```./newProject TRSM```

- copy the desired model files from DM/models into micromegas. I use h4GOn:
```cp ../models/h4GOn/* TRSM/work/models/```

-  in TRSM/ there is a main.c file and a Makefile. Compile with:
```make main=main.c``` or just ```make```
The setup is complete and and micromegas is used by codes in DM/example.

- To test manually, copy the example data point:
```cp DM/data.par micromegas_6.1.15/TRSM/```
```./main data.par```
compare the output with DM/example_test.out

## run micrOMEGAs with a steering code:
- check out DM/example_steer/README.md and create another directory for the study you'd like to perform.
