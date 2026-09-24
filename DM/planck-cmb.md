# Planck CMB constraint

New scans with `--micromegas-version 7` enable the Planck 2018 annihilation
constraint. `--no-planck-cmb` disables it. The v2 default backend is **7.1.4**, with CMB enabled. Explicit
6.1.15 runs default to CMB disabled. `--planck-cmb` explicitly enables it with either
backend, provided the TRSM driver has been rebuilt from `DM/main.c`.

## Physics prescription

The driver calls `PlanckCMB(sigmaV, SpA, SpE)` immediately after the existing
`calcSpectrum(1+2+4, ...)` call. It uses the unpropagated photon and positron
spectra. The library returns a dimensionless ratio to the 95% CL bound

\[
p_{\rm ann}<3.2\times10^{-28}\;{\rm cm^3\,s^{-1}\,GeV^{-1}}.
\]

Following the [micrOMEGAs 7 prescription, sections 5.1 and 6.3.1](https://arxiv.org/html/2606.06645v1),
the Python wrapper applies

\[
\xi=\min(1,\Omega h^2/0.12),\qquad
R_{\rm CMB}=\xi^2R_{\rm full}.
\]

`R_CMB > 1` excludes; equality passes. A valid zero annihilation signal or
zero abundance gives zero rescaled signal. Missing output, a spectrum error,
or an invalid rate, spectrum or CMB ratio is unavailable and prevents aggregate
DM acceptance. Other valid DM diagnostics remain available. The existing
relic upper cut of 0.121, SI rescaling and gamma-line rescaling are unchanged.
As with the existing wrapper, the abundance uses the driver's printed relic
density. The low-level Python API's existing `rescale=False` option instead
uses a unit abundance fraction; scan and reevaluation commands use rescaling.

Reproducible standalone commands, current v2 reference values and raw logs
are in [steer_example/README.md](steer_example/README.md).

The method is `micromegas_planck2018_swave_v1`: the built-in low-velocity,
s-wave approximation. Production settings remain spectrum key 7,
`SpectraFlag=0`, `vRot=220 km/s`, and `VZdecay=VWdecay=1` in v2.
The explicitly selected historical `TRSM_LEGACY_VIRTUAL_OFF=1` setting uses zero. The library's
default two-body annihilation velocity is `sqrt(3)*vRot/c`. This is not a
recombination-era velocity evolution calculation. Rapid velocity dependence,
resonances and final-state thresholds require care when interpreting the
approximation; the checks below do not establish its accuracy everywhere.

## Output and compatibility

Enabled scan filenames include `-cmb-planck2018`. Their eight additional
columns are `dm_cmb_enabled`, `dm_cmb_available`, `dm_cmb_status`,
`dm_cmb_reason`, `dm_cmb_ratio_raw`, `dm_cmb_abundance_fraction`,
`dm_cmb_ratio`, and `dm_cmb_excluded`. Unavailable numerical/exclusion values
are written as `nan`; an available abundance fraction can still be recorded.
The CMB status is independent of the gamma-line availability status.

Metadata and checkpoint fingerprints record the normalization, method,
backend and spectrum settings. The driver's read-only `--capabilities`
query runs before CMB-enabled scan output is created and rejects outdated
executables. Direct driver calls use `main card.dat --planck-cmb`; without
that flag its original calculations and output are retained.

Campaigns whose saved options lack CMB settings resume disabled, preserving
their columns and fingerprints, including legacy version 7 campaigns. Physics
settings cannot be changed during scan resume. Plot summaries label older
files CMB-unassessed and distinguish CMB exclusions and unavailable CMB results.

Full reevaluation is explicit and writes a separate output:

```bash
./trsmdm/bin/python reevaluate_trsm_dm_higgs.py output/old.dat \
  --output output/old_mo7_cmb.dat --micromegas-version 7
```

This recomputes physical widths, HiggsTools and all DM constraints. It
preserves the SI treatment from the input's row/metadata provenance unless
overridden with `--dm-limit-table` or `--dm-limit-model`. A missing or changed
recorded table requires an explicit replacement table. The temporary LZ 2026
curve remains explicitly approximate and is never substituted for another
selected SI treatment.

Repeat the evaluation options with `--resume`. New reevaluation checkpoints
record the backend, CMB prescription and SI table hash and reject changes.
Legacy checkpoints have no saved backend identity: resume them with their
original backend, CMB disabled and the original default SI fit. Adoption
records that configuration for subsequent resumes. Existing CMB fields are
refreshed or cleared and marked disabled. Completed inputs are never changed.

## Historical validation on manto, 2026-09-22

The results below predate the v2 SM-input, loop and virtual-W/Z corrections.
They document the original CMB implementation and are **not** the reference
values for the current example. Use its committed v2 benchmarks for reruns.

The version 7 driver was built as a separate executable before deployment.
The affected DM, SI-table, scan, resume, reevaluation, output, EWPT integration,
model-normalization and plotting tests passed in the `trsmdm` environment.
Full scan and reevaluation fixtures also exercised enabled/disabled CMB and
recovery of the selected temporary LZ 2026 SI table.

Direct calls to the installed library agree exactly with the emitted CMB
ratio in 18 benchmark outputs across ten parameter points. Selected results
use the driver's printed relic density:

| Point | MX [GeV] | Omega h² | Raw ratio | Rescaled ratio | CMB |
| --- | ---: | ---: | ---: | ---: | --- |
| Heavy | 500 | 0.266 | 0.0106474 | 0.0106474 | pass |
| Near Higgs resonance | 62 | 0.00175 | 0.0288185 | 6.129e-6 | pass |
| Light, enlarged portals | 10 | 0.119 | 1.69141 | 1.66333 | excluded |
| Higgs pole | 62.5 | 0.00325 | 11378.0 | 8.34588 | excluded |

All points have `LX=0.1` and `vevs=200 GeV`. Heavy: `Mh2=380 GeV`,
`SinT=-0.149438`, `LHX=0.05`, `LSX=0.15`. The two resonance points use
the same `Mh2` and mixing, with `LHX=0.001`, `LSX=0.01`. The enlarged
light portals are `LHX=0.28`, `LSX=0.56`, with `Mh2=80 GeV`,
`SinT=0.0998334`. These are constraint benchmarks, not all-constraints viable
points. Additional light cases bracket the CMB boundary.

Validation-only changes of `vRot` from 220 to 2.2 km/s changed the raw CMB
ratio by about 0.026% at `MX=62.5 GeV`, -0.176% at `62.49995 GeV`, and
-1.925% immediately above the W threshold at `80.3851 GeV`. At
`80.3849 GeV` the change was below 0.0002%. Production velocity and spectrum
settings were not changed. These checks quantify local sensitivity, not an
error bound or a recombination calculation.

The original and CMB-disabled version 7 drivers give identical parsed relic,
SI and gamma-line predictions at the heavy, light and near-resonance reference
points. Injection checks cover a `calcSpectrum` error, non-finite rate,
non-finite spectrum, non-finite CMB ratio, and valid zero signal. All preserve
the other parsed diagnostics; only valid zero passes the CMB check.

A separate version 6 reference installation also builds and reports a valid
CMB result (raw ratio 0.0277795 at the near-resonance point). The active
version 6 executable is unchanged. A new enabled version 7 campaign was
created and resumed using its saved settings.

Build logs, before/after drivers, benchmark cards, direct-call and fault
harnesses, raw output, `comparison.json`, regression logs and small workflow
fixtures are retained on manto under
`/Users/apapaefs/Projects/micromegas_7.1.4/validation/planck-cmb-20260922/`.
No user campaign was reevaluated or interrupted.
