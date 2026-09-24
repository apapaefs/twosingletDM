# Published direct-detection upper limits

The TRSM driver calculates elastic spin-independent (SI) nucleon cross sections
in micrOMEGAs. The Python scan applies the experimental limit afterwards, using
the neutron SI cross section in pb and the existing approximately isoscalar
Higgs-portal treatment. The same comparison works with micrOMEGAs 6.1.15 and
7.1.4; it does not change micrOMEGAs' built-in experimental likelihood routines.

The v2 default is the **LZ WS2024 observed, power-constrained 90% CL SI table**
from [HEPData 155182, version 2, table 1](https://doi.org/10.17182/hepdata.155182.v2/t1),
associated with the published [LZ analysis](https://arxiv.org/abs/2410.17036v3).
Both the production scan and `DM/steer_example/` load
[`data/lz2024/observed-si-v2.json`](data/lz2024/observed-si-v2.json) automatically.
It contains 26 mass points from 9 to 10,000 GeV and selects the `limit` column.
The original HEPData YAML, its checksum and representative numerical checks
are documented in [the table README](data/lz2024/README.md).

Verify the saved YAML checksum and every normalized table entry with:

```bash
../runtime-v2/venv/bin/python DM/tools/verify_lz2024.py
```

The inherited `lz2025-source` and `legacy-output` fits remain available through
explicit historical options in the DM provider and standalone example. New
default runs use the verified table.

## LZ September 2026: temporary approximation available

The requested paper is [arXiv:2609.02823v1](https://arxiv.org/abs/2609.02823),
"Search for dark matter particle interactions in an extended nuclear recoil
energy window with the LUX-ZEPLIN (LZ) experiment." Its linked data release is
[HEPData 182472, version 1](https://doi.org/10.17182/hepdata.182472.v1).
On 2026-09-17 that record redirected to a HEPData login page. After login, it
explicitly denied access and requested uploader, reviewer or observer
permissions. That access check concerned the September 2026 release. No official
numerical table from that release is included in this repository.

At the user's request, an optional [temporary high-mass table](data/lz2026/README.md)
uses the Figure S7 elastic upper limit at 1 TeV, approximately 4.56e-47 cm^2,
and assumes a limit proportional to DM mass from 400 to 4000 GeV. This mass
dependence is an approximation, not a digitized LZ mass curve. Select it with
`--dm-limit-table DM/data/lz2026/lz2026-figs7-highmass-approx.json`; for a random
scan also set `--m3-min 400 --m3-max 4000`. Its provenance is stored in the scan
metadata. The default remains the verified WS2024 table.

Download the complete YAML archive for the record, including `submission.yaml`
and all data tables. The needed result is the **observed upper endpoint** of the
90% confidence interval for the elastic, isoscalar scalar SI interaction as a
function of DM mass. The elastic scalar relativistic operator and its mapping
to O1 must be checked against the release's normalization. An inelastic O1
table can supply the elastic result only at delta = 0. Figure S7 alone scans
delta at a fixed mass and cannot supply a general mass-dependent scan limit.

The event mentioned in the paper remains in the published statistical analysis.
We use the upper bounds without imposing a positive signal or a lower bound;
we do not substitute the expected sensitivity or remove/reanalyse the event.

If the released quantity is Q = (c1^s * v^2)^2, supplemental Eq. (6) gives
sigma_SI^N = Q * mu_N^2 / (pi * v^4) in GeV^-2, with v = 246.2 GeV and
mu_N = m_DM * m_N / (m_DM + m_N). Conversion to pb multiplies by
3.8937966e8. The input must really be this **squared**, isoscalar O1 coefficient;
a relativistic coefficient requires its operator mapping first. The vector
conversion factor discussed below Eq. (6) is not applied to scalar TRSM DM.

## Using a verified table

The verified WS2024 table is selected automatically. To choose an alternative
release that has been checked and converted to per-nucleon cross sections,
select its normalized JSON file with:

```bash
../runtime-v2/venv/bin/python generate_trsm_points.py 123 --nrandom 500 \
  --micromegas-version 7 --dm-limit-table /path/to/verified-si-limits.json
```

The table uses schema `trsm_si_upper_limit_v1` and these required fields:

| Field | Value |
| --- | --- |
| `schema` | `"trsm_si_upper_limit_v1"` |
| `label` | Short identifier containing lowercase letters, digits, `.`, `_`, `-` |
| `source` | Publication/data DOI or URL, including the release version |
| `confidence_level` | `0.9` |
| `interaction` | `"elastic_isoscalar_si"` |
| `limit_kind` | `"observed_upper"` |
| `cross_section` | `"per_nucleon"` |
| `cross_section_unit` | `"cm2"` or `"pb"` |
| `points` | List of objects with numeric `mass_GeV` and `upper_limit` |

At least two points are required, in strictly increasing mass order, with finite
positive upper bounds. An optional `provenance` object records digitization or
approximation details and is preserved in scan metadata. A HEPData YAML export is **not** this
normalized JSON file; its columns and operator conventions must be reviewed
before conversion.

The scan converts cm^2 to pb and interpolates linearly in log(mass) and
log(cross section). Outside table coverage, DD is unassessed and the point is
retained; there is no extrapolation or fallback to another experiment.
The abundance fraction is xi = min(1, Omega/0.12): within table coverage, the SI
cross-section limit is divided by xi, with an infinite limit for xi = 0.
The separate relic-density upper cut is 0.121. Gamma-line and enabled CMB signals
use the same 0.12 abundance reference and are rescaled by xi squared.

The file's SHA-256, source, units, range and interpolation are recorded in scan
metadata; `dm_limit_model` records the label and checksum prefix. Output names
also include these identifiers. The loaded table is frozen for the run, and
resuming a campaign after changing the file is rejected. Legacy checkpoints
keep their original default-limit fingerprints.

## Historical validation on manto (2026-09-17)

These checks predate the WS2024 default and v2 abundance convention. Current
validation is documented in the [v2 readiness report](../docs/next-scan-readiness-v2.md)
and [standalone example](steer_example/README.md).

All 146 tests covering the table reader, DM checks, backend selection, scan
execution, checkpoint/resume, EWPT integration, reevaluation, output and model
normalization passed in manto's `trsmdm` environment. The original and updated
generators also gave identical fingerprints for an existing default-limit scan.
A full micrOMEGAs 7.1.4 point at MX = 62 GeV verified table selection, rescaling
and metadata with an explicitly **synthetic** table. This validates the software
connection, not the unavailable official LZ 2026 mass limits.

The source backup, test log, synthetic fixture and smoke-test manifest are under
`/Users/apapaefs/Projects/micromegas_7.1.4/validation/dd-table-support/` on manto.
