Directory layout:
`example_steer/generate_oks.py`
Generates `example_steer/oks.dat`.

`example_steer/source/`
Contains the Python scripts:
`write_mo.py`, `mO_excluder.py`

`example_steer/run/`
Contains the runtime files:
`oks.dat`, `cards/`, `MOrun.sh`, and job-submission scripts in `launch/`

`example_steer/output/`
Contains all outputs from micrOMEGAs and `mO_excluder`

`oks.dat` format:
`index LX LHX LSX MX vevs SinT Mh2`

Example scan generation:
cd DM/example_steer
python3 generate_oks.py --lhx-start 0.01 --lhx-end 1 --lhx-step 0.01 --lsx-start 0.01 --lsx-end 1 --lsx-step 0.01
It can also generate logarithmically spaced points:
python3 generate_oks.py --lhx-log --lhx-start 0.01 --lhx-end 1 --lhx-num-points 100 --mx-start 10 --mx-end 1000 --mx-step 10
To set a default value for a parameter x and not vary it, use --x-start and do not provide --x-end.

Constrained scan generation:
`--equal-couplings` writes `LSX = LHX` for every generated point. The `LHX` scan options define the shared coupling values, and the `LSX` scan options are ignored while this rule is active.

Example:
python3 generate_oks.py --equal-couplings --lhx-start 0.01 --lhx-end 1 --lhx-step 0.01

`--resonant-mx` writes `MX = Mh2 / 2` for every generated point. The `Mh2` scan options define the resonance values, and the `MX` scan options are ignored while this rule is active.

Example:
python3 generate_oks.py --resonant-mx --mh2-start 200 --mh2-end 1000 --mh2-step 100

The two rules can be combined:
python3 generate_oks.py --equal-couplings --resonant-mx --lhx-start 0.01 --lhx-end 1 --lhx-step 0.01 --mh2-start 200 --mh2-end 1000 --mh2-step 100

Copy `oks.dat` into the run directory:
cp oks.dat run/oks.dat

Run the Python scripts (Python 3.6+ required):
To generate micrOMEGAs cards:
cd source
python3 write_mo.py

To see the detection limits used, plot them using:
python3 mO_excluder.py --plot-dirdet-limits
python3 mO_excluder.py --plot-indirect-limits

Run the full micrOMEGAs + exclusion chain:
cd ../run
MICROMEGAS_MAIN=<full_path>/twosingletDM/DM/micromegas_6.1.15/TRSM/main ./MOrun.sh
(the default MICROMEGAS_MAIN in MOrun.sh should also work if you run ./MOrun.sh)

What the workflow writes:
`run/oks.dat`
Copy of the scan definition used for the run.

`run/cards/MO_inp*.dat`
Input cards written by `source/write_mo`.

`output/OUT_mO_<index>`
Raw micrOMEGAs output for each processed point.

`output/MOrun.log`
Always generated. Captures the terminal output printed by `run/MOrun.sh`.

`output/scan_results.dat`
Always generated. One line per processed point:
`index LX LHX LSX MX vevs SinT Mh2 MDM Omega DirDet`

`output/DM_data/DM_data_<index>`
Always generated per processed point. Summary from `source/mO_excluder`:
`index LX LHX LSX MX vevs SinT Mh2 MDM Omega DirDet DirDetLimit DirDetBaseLimit IndirAvailable IndirEnergy IndirFlux IndirLimit IndirRatio`

`output/allall.dat`
Generated anyway but it is not empty if at least one point passes the relic-density upper limit,
direct-detection limit, and indirect-detection limit.
Accepted points with the same columns as `DM_data_<index>`.

`output/all_dirpass.dat`
Generated anyway but it is not empty if at least one point passes the relic-density upper limit and
the direct-detection limit.
This includes points that later fail indirect-detection checks.
Columns are the same as `DM_data_<index>`.

`output/relic_pass.dat`
Generated anyway but it is not empty if at least one point passes the relic-density upper limit
`Omega <= 0.121`. Columns are the same as `DM_data_<index>`.

`output/relic_strict.dat`
Generated anyway but it is not empty if at least one point satisfies the strict relic-density band
`0.119 <= Omega <= 0.121`. Columns are the same as `DM_data_<index>`.

`output/omexcl.dat`
Generated anyway but it is not empty if at least one point is excluded by relic density.
Columns:
`index LX LHX LSX MX vevs SinT Mh2 MDM Omega DirDet DirDetLimit DirDetBaseLimit IndirAvailable IndirEnergy IndirFlux IndirLimit IndirRatio`

`output/luxexcl.dat`
Generated anyway but it is not empty if at least one point is excluded by direct detection.
Columns:
`index LX LHX LSX MX vevs SinT Mh2 MDM Omega DirDet DirDetLimit DirDetBaseLimit IndirAvailable IndirEnergy IndirFlux IndirLimit IndirRatio`
and `output/luxexcl.dat` is generated otherwise.

points that pass Relic Density limit but fail direct detection are in `omgpass_dirfail.dat`


`output/indirexcl.dat`
Generated anyway but it is not empty if at least one point is excluded by indirect detection.
The indirect check reads `FermiLAT_line_channel` lines from the current
micrOMEGAs output, compares each channel's `Phi_R16` integrated line flux
against the interpolated Fermi-LAT R16 gamma-line flux limit from
arXiv:1506.00013, and excludes the point if any channel has
`IndirRatio = IndirFlux / IndirLimit > 1`.
Columns:
`index LX LHX LSX MX vevs SinT Mh2 MDM Omega DirDet DirDetLimit DirDetBaseLimit IndirAvailable IndirEnergy IndirFlux IndirLimit IndirRatio`

`output/indirpass.dat`
Generated anyway but it is not empty if at least one point is not excluded by indirect detection.
so `output/indirpass.dat` + `output/indirexcl.dat` = all points.
This means that some (most) of the points in `output/indirpass.dat` are excluded by other limits.
Columns:
`index LX LHX LSX MX vevs SinT Mh2 MDM Omega DirDet DirDetLimit DirDetBaseLimit IndirAvailable IndirEnergy IndirFlux IndirLimit IndirRatio`

If a point passes relic density limit and direct detection but does not pass indirect detection,
it is detected by indirect detection so it is saved in `indir_caughtit.dat`.
In the same manner, there is `dir_caughtit.dat` and `dir_indir_caughtit.dat` when a point passes
relic density limit but is detected both directly and indirectly.

`output/dmexcl.dat`
Generated anyway but it is not empty if at least one point fails any DM check.
so allall.dat + dmexcl.dat = scan_results.dat
Columns are the same as `DM_data_<index>`.

Marker files (commented out in mO_excluder for now):
`output/DM_EXCLUDED/DM_EXCLUDED_<index>`, `output/RelDens_EXCLUDED/RelDens_EXCLUDED_<index>`,
`output/DirDet_EXCLUDED/DirDet_EXCLUDED_<index>`, `output/IndirDet_EXCLUDED/IndirDet_EXCLUDED_<index>`
These are created only for points that fail the corresponding check.

Plotting utilities:
`plot/plot_omega.py`
Reads `output/scan_results.dat` and plots `Omega` against one chosen variable.

Example:
python3 plot/plot_omega.py output MX

This writes `plot/outplots/omega_vs_MX.png` and `plot/outplots/log10_omega_vs_MX.png` by default.

`plot/plot_omega_colored.py`
Reads `output/scan_results.dat` and plots `Omega` against one chosen variable,
with another variable given a color map.

Example:
python plot_omega_colored.py /path/to/output MX LSX

This writes `plot/outplots/omega_vs_MX_colored_by_LSX.png` and `plot/outplots/log10_omega_vs_MX_colored_by_LSX.png` by default.

`plot/steer_plots.py`
Reads `run/oks.dat`, detects which scan variables vary, and runs the plotting
scripts automatically.

Default usage:
`python3 plot/steer_plots.py`

This uses `output/` automatically as the output directory. If `MX` is one of the
varying variables in `run/oks.dat`, it is used as the x-axis. The second varying
variable is used as the y-axis for the 2D relic plots and as the color variable
for `plot_omega_colored.py`.

The driver also detects constrained scans in `run/oks.dat`. If both `MX` and
`Mh2` vary while all rows satisfy `MX = Mh2 / 2`, it keeps `MX` on the x-axis,
avoids using `Mh2` as the secondary variable when another varying variable is
available, and appends `(MX = Mh2 / 2)` to plot titles. If both `LHX` and `LSX`
vary while all rows satisfy `LHX = LSX`, it uses `LHX` as the secondary variable
when possible and appends `(LHX = LSX)` to plot titles.

If only one scan variable varies, the driver still runs `plot_omega.py` and skips
the plots that require a second variable.

Example with explicit overrides:
`python3 plot/steer_plots.py --outdir output_1 --xvar MX --yvar LSX`

Options:
`--outdir <directory>` overrides the default output directory.
`--xvar <variable>` overrides the automatically selected x-axis variable.
`--yvar <variable>` overrides the automatically selected secondary variable.


`plot/plot_relic_pass_2d.py`
Reads `output/scan_results.dat`, keeps only points with `Omega <= 0.1224`, and
plots two chosen variables against each other.

Example:
python3 plot/plot_relic_pass_2d.py output MX LSX

This writes `plot/outplots/relic_pass_LSX_vs_MX.png` and `plot/outplots/relic_pass_log10_LSX_vs_MX.png` by default.

`plot/plot_relic_strict_2d.py`
Reads `output/scan_results.dat`, keeps only points with
`0.121 <= Omega <= 0.119`, and plots two chosen variables
against each other.

Example:
python3 plot/plot_relic_strict_2d.py output MX LSX

This writes `plot/outplots/relic_strict_LSX_vs_MX.png` and `plot/outplots/relic_strict_log10_LSX_vs_MX.png` by default.

All plotting scripts accept `--output` to choose the image filename and `--show`
to display the plot interactively.
