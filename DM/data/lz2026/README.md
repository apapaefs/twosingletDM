# Temporary LZ 2026 high-mass limit

`lz2026-figs7-highmass-approx.json` is an **approximation**, not the official
HEPData mass table. It is available as an explicit scan option while that
release is inaccessible. The existing default fit is unchanged.

The input is the solid observed upper curve of Figure S7 in
[arXiv:2609.02823v1](https://arxiv.org/abs/2609.02823v1), at delta = 0 and a
DM mass of 1,000 GeV. Reading its vector path gives 4.5555509e-47 cm^2, rounded
to **4.56e-47 cm^2** for the temporary table. The published interval includes
the event discussed in the paper; no event is removed and no positive-signal
or lower-limit condition is imposed.

The selected approximate mass dependence is

```text
sigma_limit(m) = 4.56e-47 * (m / 1000 GeV) cm^2,  400 <= m/GeV <= 4000.
```

| DM mass [GeV] | Approximate SI upper limit [cm^2] | [pb] |
| ---: | ---: | ---: |
| 400 | 1.824e-47 | 1.824e-11 |
| 1000 | 4.56e-47 | 4.56e-11 |
| 4000 | 1.824e-46 | 1.824e-10 |

This uses the heavy-WIMP, approximately fixed spectral shape scaling
rate ~ sigma/m. The paper's look-elsewhere discussion notes nearly degenerate
recoil spectra for masses above 400 GeV, but **does not establish this formula
for the limits**. Changes in recoil shape and likelihood response are neglected.
No uncertainty or confidence coverage has been established for the scaling
away from the plotted 1 TeV anchor. In particular, this approximation must not
be extended to the low-mass or Higgs-resonance regions.

From the scan repository on manto:

```bash
./trsmdm/bin/python generate_trsm_points.py 123 --nrandom 500 \
  --micromegas-version 7 \
  --dm-limit-table DM/data/lz2026/lz2026-figs7-highmass-approx.json \
  --m3-min 400 --m3-max 4000
```

The same table can be selected with version 6. The existing relic-density
rescaling is applied after this base limit. The scan checks mass coverage
before writing output, records the approximation and source in its metadata,
and includes `highmass-approx` and a checksum in the run name. Omitting
`--dm-limit-table` restores the existing `lz2025-source` treatment.

## Reproducing the extraction

Download and unpack the [arXiv v1 source archive](https://arxiv.org/src/2609.02823v1),
then run (requires `pdfplumber`):

```bash
python DM/tools/extract_lz2026_figs7.py /path/to/extracted-arxiv-source --check
```

Omit `--check` to regenerate the two JSON files. The extractor verifies the
SHA-256 of the original figure PDFs and selects the eight vertices of the
solid black **upper** curve. It excludes the dashed median sensitivity and
the solid lower bound. Figure S7's horizontal axis is delta, not DM mass.
`figure-s7-o1-1tev-digitized.json` preserves all eight extracted splitting
points and their PDF coordinates as an audit record; it uses a different
schema and cannot be passed to `--dm-limit-table`.

As a normalization cross-check, the script independently extracts the squared
O1 coefficient from Figure 6's top panel and converts it with supplemental
Eq. (6), v = 246.2 GeV and m_N = 0.939 GeV. The cross sections agree with
Figure S7 within **0.467%**. The table uses Figure S7 directly. This agreement
checks the extraction/normalization; it is not an uncertainty on the assumed
mass dependence. Replace the approximation with verified HEPData values when
they become available.

## Validation on manto, 2026-09-17

All 151 regression tests passed. An actual micrOMEGAs 7.1.4 scan point at
M3 = 500 GeV used the interpolated base limit of 2.28e-11 pb and recorded the
approximation provenance. That benchmark was excluded by its relic density
and SI cross section, as expected. A zero-draw campaign resumed with the same
table, and a range starting at 399 GeV was rejected before creating output.
The extraction also reproduced both JSON files exactly from the source PDFs.

Source backups and validation logs are retained on manto under
`/Users/apapaefs/Projects/micromegas_7.1.4/validation/lz2026-temporary/`.
