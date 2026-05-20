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
- To begin the scan, execute ```generate_trsm_points.py SEED``` where the SEED is an integer to be used as a seed for the random numbers. 

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
