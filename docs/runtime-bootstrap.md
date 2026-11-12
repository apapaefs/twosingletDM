# Automated runtime installation

`tools/bootstrap_runtime.py` downloads and builds the software used by the
active TRSM pipeline in a separate directory. Use the same commands on the
laptop, Manto, or a Linux host. Run from the repository checkout with Python
3.13; no machine-specific paths are embedded in the installer.

Preview the downloads and builds without changing any files:

```bash
python3.13 tools/bootstrap_runtime.py --prefix ../runtime-v2-new --dry-run
```

Build the complete runtime and activate it:

```bash
python3.13 tools/bootstrap_runtime.py --prefix ../runtime-v2-new --jobs 4
source ../runtime-v2-new/activate.sh
```

The first invocation requires a fresh directory. An existing unmanaged
installation is refused, including the current `runtime-v2`. Activation selects
the new Python environment, native executables, Higgs datasets, and generated
MadGraph processes. It does not modify shell startup files or launch a scan.
The old runtime remains available for existing campaigns.

## What is installed

| Component | Source and preparation |
| --- | --- |
| Python packages | An isolated virtual environment with `config/runtime-requirements-v2.txt` and the build tools in `config/build-requirements-v2.txt` |
| HiggsTools 1.1.3 | The recorded upstream revision, including pinned CMake FetchContent dependencies and compatibility flags for current macOS compilers |
| HiggsBounds and HiggsSignals data | The dataset revisions used by the validated v2 runtime, cloned inside the new prefix |
| micrOMEGAs 7.1.4 and 6.1.15 | Checksum-verified release archives, finite-value patches, canonical TRSM model/loop code, and compiled `TRSM/main`; each build is serial because upstream static libraries are shared |
| BSMPT 3.2.1 | Public upstream revision `431df3b6…`, the repository's TRSM model and precision helpers, shared SM inputs, Conan dependencies, and `CalcTemps`, `MinimaTracer`, `PhaseProbe`, `Test` |
| MadGraph 3.5.15 | Pinned upstream revision, the tracked `MG5stuff/loop_sm_twoscalar_generic.tar.gz` UFO export, COLLIER 1.2.9 and bundled CutTools/IREGI, and generated/compiled process directories |

The five generated processes match `generate_mg5_trsm_xsecs.ProcLocation`:

| Scan key | MadGraph process | Directory |
| --- | --- | --- |
| `hh` | `g g > h h [noborn=QCD]` | `gg_hh_twoscalar` |
| `hhh` | `g g > h h h [noborn=QCD]` | `gg_hhh_twoscalar` |
| `gg_heta0` | `g g > h eta0 [noborn=QCD]` | `gg_heta0` |
| `pp_eta0Z` | `p p > eta0 z` | `pp_eta0Z` |
| `gg_eta0Z` | `g g > eta0 z [noborn=QCD]` | `gg_eta0Z` |

Generation and compilation produce no event samples. Process cards retain the
bundled `nn23lo1` central PDF prescription. The optional LHAPDF uncertainty pass
is disabled; installing LHAPDF and selecting an uncertainty prescription is a
separate analysis choice. COLLIER and the bundled reduction libraries suffice
for these processes; optional Ninja, showers, detector simulation, and analysis
packages are not installed. CalcHEP is included with micrOMEGAs. The canonical
generated CalcHEP tables and UFO are already tracked, so a build does not require
LanHEP or Mathematica/FeynRules or regenerate physics model exports.
The initial tree-process executable uses the original matrix elements;
MadGraph can prepare its usual helicity optimization when a run is launched.
For loop-induced processes the installer also compiles each MadLoop `check`
executable and verifies its local Fortran modules. Building the integration
executable alone is insufficient: the parent-directory `OLP_static` build can
leave `polynomial.o` inside `PV*/` but its `.mod` file in `SubProcesses/`.
The installer forces a local polynomial compilation when that module is missing,
then records both the module and the `check` executable as installation artifacts.
Fresh exports also retain MG5's reduction-library runtime paths in the MadLoop
check link command, so a shared COLLIER library can be loaded at initialization.

## Adding a process to an existing runtime

After activating the existing runtime, add the loop-induced ZH2 process with:

```bash
python tools/setup_mg5_process.py --mg5-process gg_eta0Z
```

The command uses `TRSM_MG5_LOCATION` (or `--mg5-location PATH`), the installed
UFO and COLLIER, and the same central-PDF compilation procedure as bootstrap.
It records a process/model receipt and build logs, validates an existing
matching process, and refuses to replace a conflicting process. It does not
save changes to the global MG5 configuration. Bootstrap `--resume` cannot add
a process because its original process list is part of the installation identity.
See the [saved-campaign commands](../README.md#add-loop-induced-zh2-rates-to-a-saved-campaign)
to augment existing points separately from their original scan checkpoints.

## Repairing an existing MadLoop installation

If a previous installation fails with `Cannot open module file
'ml5_..._polynomial_constants.mod'`, keep its process directories and saved
campaigns. From the repository root, after activating that runtime, check the
selected processes without modifying them:

```bash
python tools/madloop_build.py \
  "$TRSM_MG5_LOCATION/gg_heta0" \
  "$TRSM_MG5_LOCATION/gg_hh_twoscalar" \
  "$TRSM_MG5_LOCATION/gg_hhh_twoscalar"
```

With scans using those processes stopped, rebuild the missing artifacts:

```bash
python tools/madloop_build.py --repair \
  "$TRSM_MG5_LOCATION/gg_heta0" \
  "$TRSM_MG5_LOCATION/gg_hh_twoscalar" \
  "$TRSM_MG5_LOCATION/gg_hhh_twoscalar"
```

Select only directories that exist in your installation. This command uses the
same exclusive process lock as the scan interface and refuses an active writer.
It runs serial `make` builds and keeps a new `madloop-build-*.log` for every
repair attempt. It preserves source files, cards, saved events, and campaign
checkpoints, and does not run `check`, generate events, or restart a campaign.
Tree-level processes can be checked too; they need no MadLoop repair.

After repair, resume a campaign with its original checkout, environment, and
settings. The repair changes only compilation artifacts and therefore preserves
its runtime receipt. Updating the scan's Python sources or Git revision can
invalidate its fingerprint; do not edit the saved fingerprint to bypass this
protection. The repair script uses only the Python standard library and can be
run from a separate checkout while retaining the original campaign checkout.

## Repairing an existing MG5 survey

Older scan wrappers requested one survey iteration for loop-induced processes,
but MG5 3.5.15's unsplit integration requires at least three. This can produce
`Survey return zero cross section` even when individual integration channels
have nonzero estimates. With three iterations, its adaptive helicity sampler
can also fail with `DiscreteSampler:: Error, no point could be picked`.
New wrappers request three iterations and explicitly sum loop helicities
(`nhel=0`). This preserves the unpolarized process and physics parameters.

To recover saved campaigns using their original checkout, run this maintenance
tool from the updated checkout while all scans using these processes are stopped:

```bash
# Check first; --repair applies the changes.
python tools/repair_mg5_survey.py \
  "$TRSM_MG5_LOCATION/gg_heta0" \
  "$TRSM_MG5_LOCATION/gg_hh_twoscalar" \
  "$TRSM_MG5_LOCATION/gg_hhh_twoscalar"
python tools/repair_mg5_survey.py --repair \
  "$TRSM_MG5_LOCATION/gg_heta0" \
  "$TRSM_MG5_LOCATION/gg_hh_twoscalar" \
  "$TRSM_MG5_LOCATION/gg_hhh_twoscalar"
```

Select only installed directories. The tool uses the scan's process writer lock
and refuses an active writer. In each generated process it guards the survey
input writer against `maxiter < miniter`, logging when it raises that invalid
request. Legitimate split-grid jobs with `maxiter=miniter=1` are preserved.
For loop-induced processes it also sets `nhel=0` in the current and default run
cards. All other card settings, parameter cards, matrix elements, UFO sources,
events, and campaign files are preserved. Every changed file has a backup and
SHA-256 receipt under the process's `.trsm-maintenance/survey-minimum-*/` directory.
Repeated repair is a no-op unless the settings need repair again.

This repair is deliberately limited to the generated runtime; it does not edit
the original scan checkout or rewrite saved fingerprints. Existing receipts
exclude these mutable integration files, so the original campaign can resume
with its original checkout and options. Keep the maintenance receipts alongside
the campaign's provenance. A future process regeneration needs the repair again
if it is used with an old wrapper. The current wrapper supplies the corrected
settings on every invocation. The repair command does not restart scans.

## Prerequisites

Install Python 3.13, Git, C/C++ and Fortran compilers, `make`, `patch`, and
`rsync`. CMake, Conan, and the Ninja build tool are installed into the virtual
environment. The Python Ninja build tool is distinct from the optional
MadLoop reduction library of the same name.

On macOS, install the Xcode command-line tools and use Homebrew, for example:

```bash
xcode-select --install
brew install python@3.13 gcc git rsync
```

On Debian/Ubuntu, install `build-essential`, `gfortran`, `git`, `patch`,
`rsync`, and `pkg-config` through the system package manager, and provide
Python 3.13 with its `venv` and development headers. The installer reports
missing prerequisites before downloading. Compiler executables may be selected
with `CC`, `CXX`, and `FC`. Native build paths must not contain whitespace.
System package installation is left to the host administrator; the build itself
does not use `sudo`.

## Restarting and selecting components

To continue after a network or compiler failure, repeat the original command
with `--resume`:

```bash
python3.13 tools/bootstrap_runtime.py --prefix ../runtime-v2-new --jobs 4 --resume
```

Completed artifacts are verified before a step is skipped. A failed step's
partial directories are moved under `failed-builds` before retrying. Logs are
appended under `logs`, and only one installer may use a prefix at a time.
Changed source pins, model inputs, compiler selections, component selection, or installer
implementation require a new prefix. `--jobs` may change on resume.

The first loop-process export can spend several minutes compiling CutTools and
IREGI without printing another progress line. Each step prints its log path;
use `tail -f` on that file to inspect it. `--jobs` controls parallel CMake/Conan
builds and MadGraph's core setting. The micrOMEGAs and generated-process
Makefiles run serially to avoid shared-library races.

For a DM-only installation, or a smaller MadGraph installation:

```bash
python3.13 tools/bootstrap_runtime.py --prefix ../runtime-dm-new --components micromegas
python3.13 tools/bootstrap_runtime.py --prefix ../runtime-mg5-new --components mg5 \
  --mg5-process gg_heta0 --mg5-process pp_eta0Z
```

`--components` also accepts `python`, `datasets`, `higgstools`, and `bsmpt`;
HiggsTools automatically includes its datasets. `--download-only` fetches the
selected sources without compiling. Continue with the same selection and
`--resume`, omitting `--download-only`. `--download-cache PATH` shares verified
archives across fresh installations; a corrupt cached archive is rejected.

## Provenance and checks

Source commits, release URLs, archive checksums, and process definitions live
in `config/runtime-sources-v2.json`. They come from the official
[HiggsTools repositories](https://gitlab.com/higgsbounds),
[BSMPT repository](https://github.com/BSMPT/BSMPT),
[MadGraph repository](https://github.com/mg5amcnlo/mg5amcnlo),
[micrOMEGAs releases](https://micromegasdm.github.io/), the authors'
[archived 6.1.15 release](https://zenodo.org/records/13376690), and
[COLLIER downloads](https://collier.hepforge.org/downloads/).
BSMPT's upstream Conan recipe specifies its native dependency versions; build
tool versions and Conan resolutions are recorded in the installed files/logs.
The installer does not claim byte-identical binaries across compilers or hosts.

The final checks include Python dependency consistency, flavour-data validation,
loading HiggsTools and both datasets, micrOMEGAs capability probes, BSMPT CLI
startup, and existence of compiled MadGraph subprocesses and local MadLoop
initialization artifacts. These are installation
checks, not a new physics-validation campaign.
BSMPT 3.2.1's help command returns status 1 after its missing-argument message;
the startup check accepts that specific response only when the expected help
and the registered TRSM model are present.

`setup-manifest.json`, `bootstrap-state.json`, step logs, and `pip-freeze.txt`
record the inputs and artifacts. A complete installation additionally runs
`tools/runtime_manifest.py`, which verifies native sources against this checkout
before writing the pipeline's `runtime-manifest.json`. Higgs dataset paths are
selected with `TRSM_HB_DATASET` and `TRSM_HS_DATASET`; the provider and physics
manifest use the same paths. Without these variables, historical sibling
dataset locations remain the defaults.
