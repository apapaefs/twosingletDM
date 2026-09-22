# micrOMEGAs 7.1.4 on manto

Installed and validated on 2026-09-17. The existing 6.1.15 installation and
running reevaluation were left in place.

## Installation

- New executable: `/Users/apapaefs/Projects/micromegas_7.1.4/TRSM/main`.
- Existing executable: `/Users/apapaefs/Projects/micromegas_6.1.15/TRSM/main`.
- Official archive: <https://micromegasdm.github.io/downloadarea/v7.1/micromegas_7.1.4.tgz>.
- SHA-256: `c8cf207b17541a5b7d7e7ff157f1c1eb36e1dd3f7547e7745559b195c2a34264`.
- Built with `sh DM/setup_micromegas.sh ../micromegas_7.1.4` and manto's
  existing clang/gfortran configuration. The executable is native arm64.
- `DM/main.c` and all five `DM/models/h4GOn/*.mdl` files match the installed
  6.1.15 files byte for byte. The finite-value `calcSpectrum` guard was ported
  to 7.1.4; no other upstream source changes were made.

## Validation

The 117 tests covering version selection, scan execution, checkpoint/resume,
DM parsing and constraints, model normalization, reevaluation, and scan output
passed in manto's `trsmdm` environment. The previous generator and the updated
generator produced identical configuration fingerprints for a representative
legacy scan using the default 6.1.15 backend.

A separate 6.1.15 installation under the new installation's `validation/reference`
directory supplied the reference calculations. Its C/header files in `sources`,
`include`, and `CalcHEP_src/c_source` match the active 6.1.15 installation
(excluding the installation-path header). This avoids concurrent use of the
active campaign's CalcHEP workspace.

The following card inputs were used, with `SinT=sin(a12)` and `LX=0.1`:

| Point | M2 [GeV] | MX [GeV] | vs [GeV] | a12 | LHX | LSX |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Heavy | 380 | 500 | 200 | -0.15 | 0.05 | 0.15 |
| Light | 80 | 10 | 200 | 0.1 | 0.01 | 0.02 |
| Near Higgs resonance | 380 | 62 | 200 | -0.15 | 0.001 | 0.01 |

Results below use the driver's printed precision:

| Point | Omega h², 6.1.15 | Omega h², 7.1.4 | Neutron SI [pb], both | DM constraints, both |
| --- | ---: | ---: | ---: | --- |
| Heavy | 0.266 | 0.266 | 1.526e-10 | Fail |
| Light | 63.4 | 58.6 | 1.174e-8 | Fail |
| Near Higgs resonance | 0.00175 | 0.00175 | 9.667e-12 | Pass |

The light-point relic density decreases by about 7.6% in 7.1.4. The R16
gamma-line fluxes agree to about 6 parts per million at the two points within
the line-limit energy range. These three checks establish working integration,
not agreement throughout the parameter space.

The near-resonance point also passed every enabled check in a full explicit
`generate_trsm_points.py` run with `--micromegas-version 7`. Its scan filename
contains `-mo7.1.4`, and its metadata records the version and absolute executable.
A zero-draw version 7 campaign was created and resumed without repeating any
backend options.

Logs, raw benchmark outputs, `comparison.json`, build settings, and backups of
the previous scan sources are retained on manto under
`/Users/apapaefs/Projects/micromegas_7.1.4/validation/`. The build log is
`/Users/apapaefs/Projects/micromegas_7.1.4-trsm-build.log`.

An initial smoke test using the existing `--write-all-points` mode also tried
to run EWPT and reported a missing configured BSMPT `MinimaTracer` executable.
This was a hardcoded laptop path; both binaries were installed on manto.
The follow-up fix on 2026-09-17 resolves the sibling BSMPT build relative to
the scan repository. It uses
`/Users/apapaefs/Projects/BSMPT/build/macos-armv8-release/bin` on manto and
preserves the laptop's existing build location.

The original heavy point was rerun with `--write-all-points` and no executable
overrides. Both BSMPT programs completed, and the scan recorded
`ewpt_status=success` with an empty `ewpt_error`. The 121 EWPT, scan,
reprocessing, checkpoint, and backend-selection regression tests also passed.
The small CalcTemps parser fixture is now included in the repository, so these
tests do not depend on a generated output file in the external BSMPT tree.
The follow-up logs and pre-edit source backup are under `validation/ewpt-path-fix/`.
