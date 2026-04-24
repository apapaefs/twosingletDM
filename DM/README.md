Directory layout:
`dmonly/generate_oks.py`
Generates `dmonly/oks.dat`.

`dmonly/source/`
Contains the C++ sources and compiled executables:
`write_mo.cpp`, `write_mo`, `mO_excluder.cpp`, `mO_excluder`

`dmonly/run/`
Contains the runtime files:
`oks.dat`, `cards/`, `MOrun.sh`, and job-submission scripts in `launch/`

`dmonly/output/`
Contains all outputs from micrOMEGAs and `mO_excluder`

`oks.dat` format:
`index LX LHX LSX MX vevs SinT Mh2`

Example scan generation:
cd /project/atlas/users/okarkout/theo/dmonly
python3 generate_oks.py --mx-start 10 --mx-end 1000 --mx-step 10 --lsx-start 0.01 --lsx-end 1 --lsx-step 0.01

Copy `oks.dat` into the run directory:
cp oks.dat run/oks.dat

Compile the C++ tools in `source/`:
cd source
g++ -O2 -std=c++11 -o write_mo write_mo.cpp
g++ -O2 -std=c++11 -o mO_excluder mO_excluder.cpp

Generate micrOMEGAs cards in `run/cards/`:
./write_mo

Run the full micrOMEGAs + exclusion chain:
cd ../run
MICROMEGAS_MAIN=/project/atlas/users/okarkout/theo/micromegas_6.1.15/TRSM/main ./MOrun.sh
(./MOrun.sh should also work)

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

`output/DM_data_<index>`
Always generated per processed point. Summary from `source/mO_excluder`:
`index LX LHX LSX MX vevs SinT Mh2 MDM Omega DirDet DirDetLimit LUXBaseLimit`

`output/allall.dat`
Generated only if at least one point passes all DM checks.
Accepted points with the same columns as `DM_data_<index>`.

`output/omexcl.dat`
Generated only if at least one point is excluded by relic density.
Columns:
`index LX LHX LSX MX vevs SinT Mh2 MDM Omega`

`output/luxexcl.dat`
Generated only if at least one point is excluded by direct detection.
Columns:
`index LX LHX LSX MX vevs SinT Mh2 MDM Omega DirDet DirDetLimit LUXBaseLimit`

`output/dmexcl.dat`
Generated only if at least one point fails any DM check.
Columns:
`index LX LHX LSX MX vevs SinT Mh2 MDM Omega DirDet DirDetLimit LUXBaseLimit`

Marker files:
`output/DM_EXCLUDED_<index>`, `output/RelDens_EXCLUDED_<index>`, `output/DirDet_EXCLUDED_<index>`
These are created only for points that fail the corresponding check.

Plotting utilities:
`plot/plot_omega.py`
Reads `output/scan_results.dat` and plots `Omega` against one chosen variable.

Example:
python3 plot/plot_omega.py output MX

This writes `plot/outplots/omega_vs_MX.png` and `plot/outplots/log10_omega_vs_MX.png` by default.

`plot/plot_relic_pass_2d.py`
Reads `output/scan_results.dat`, keeps only points with `Omega <= 0.1224`, and
plots two chosen variables against each other.
Produces a heatmap of Omega values and highlights in purple points that pass
the strict relic abundance criterium `0.1199 - 0.025 <= Omega <= 0.1199 + 0.025`

Example:
python3 plot/plot_relic_pass_2d.py output MX LSX

This writes `plot/outplots/relic_pass_LSX_vs_MX.png` and `plot/outplots/relic_pass_log10_LSX_vs_MX.png` by default.

`plot/plot_relic_strict_2d.py`
Reads `output/scan_results.dat`, keeps only points with
`0.1199 - 0.025 <= Omega <= 0.1199 + 0.025`, and plots two chosen variables
against each other.

Example:
python3 plot/plot_relic_strict_2d.py output MX LSX

This writes `plot/outplots/relic_strict_LSX_vs_MX.png` and `plot/outplots/relic_strict_log10_LSX_vs_MX.png` by default.

All plotting scripts accept `--output` to choose the image filename and `--show`
to display the plot interactively.
