# Upsilon lepton-product bounds

New scans apply `upsilon_leptons_v1`, using the collaborator's
`flavour/CheckUpsLeptonBounds` prescription for
`BR(Upsilon(1S) -> gamma hi) * BR(hi -> mu+mu- / tau+tau-)`.
The data are digitizations of Belle and BaBar curves in Fig. 4(c,d) of
[arXiv:2112.11852](https://arxiv.org/pdf/2112.11852).
The production correction, physical constants, and four CSV files are retained.
These are prompt, narrow-resonance bounds; no additional lifetime or broad-width
acceptance cut is imposed. The stronger available individual bound is applied
for each channel. This is not a combined experimental likelihood.

## Decays and coverage

The doublet component `ki = R[i,0]` is passed to the checker, which applies
`abs(ki)**2`. Base BR array entries 2 and 1 are respectively muons and taus.
Each is multiplied by `base_width / physical_total_width`, including exotic
and invisible partial widths. HiggsTools' provider width floor is never used.
The underlying base arrays are not modified.

The search windows are inclusive `[0.212, 9.2]` GeV for muons and
`[3.55, 9.2]` GeV for taus. Each experiment additionally requires CSV coverage.
The supplied tau curves actually start at about 3.834 GeV (Belle) and
3.839 GeV (BaBar). Interpolation is linear in mass and logarithmic in the limit;
there is no extrapolation. Equality at an upper limit passes.

The existing SM-like decay provider starts at 4 GeV, so nonzero-mixing scalars
below 4 GeV inside the experimental window remain unassessed. There is no new
low-mass hadronic decay model. Zero doublet projection has zero visible signal
without needing a decay table. In the active `vx=0` branch, `h3` is stable and
has zero projection; `h1` lies above coverage. The relevant decaying scalar is
`h2`, with possible dilution through `h2 -> h3 h3`.

## Output and selections

The additive v2 output columns are:

| Column | Meaning |
| --- | --- |
| `flavour` | `True`, `False`, or `nan` for an unassessed result |
| `flavour_method` | `upsilon_leptons_v1` |
| `flavour_status` | `passed`, `excluded`, `outside_coverage`, `zero_signal`, or `unassessed` |
| `flavour_reason` | Scalar/channel explanation, including missing-input errors |
| `flavour_max_ratio` | Largest covered prediction/limit ratio, or `nan` when none applies |
| `flavour_details` | JSON keyed by scalar, containing decay inputs and channel diagnostics |

Outside coverage is a non-veto (`True`), not experimental confirmation. A known
exclusion wins over an unavailable scalar assessment. Otherwise an unavailable
required assessment prevents a pass. JSON diagnostics use `null`, not nonstandard
NaN/Infinity, for unavailable numbers and bounds.

Full viability now requires tree theory, the existing `experimental_subset`,
`flavour is True`, and `dm is True`. Non-DM viability and `--mg5-without-dm`
drop only the DM requirement. Flavour does not change `experimental_subset`,
`ewpt_eligible`, or EWPT candidate definitions. Every evaluated scan row is
retained. Campaign counts include the independent flavour verdict.

## Saved scans and reproducibility

Run from the repository using its Python environment:

```bash
python reevaluate_trsm_flavour.py INPUT --output OUTPUT
python plot_trsm_constraint_suite.py OUTPUT --format both
```

The flavour-only command supports `vx=0` TSVs, including `.dat` files with TSV
contents. It reconstructs the relevant decays from stored `M2`, `M3`, `a12`,
`vs`, `lPhiX`, and `lSX`; `k2` is derived from the angle, so blank cached mixing
and obsolete widths do not enter the calculation. Within the Upsilon window
`h2 -> h1 h1` is closed. The existing SM decay and canonical portal helpers
supply all possible `h2` widths in this branch.

Only flavour columns are added/replaced. Other verdicts, widths, and fields
remain historical; they have not been revalidated. Unsupported/missing point
inputs produce unassessed flavour rows. Malformed TSV structure aborts before
publishing the output. The command refuses in-place writes and existing output
paths, verifies the input did not change, and records its hash, current flavour
inputs, status counts, and inherited metadata in a separate sidecar.

Generation and ordinary DM/Higgs reevaluation also calculate the same flavour
diagnostics. Source files, constants, CSVs and decay inputs are fingerprinted.
Pre-flavour or otherwise incompatible scan checkpoints are rejected before
appending. A flavour-only output is a derived dataset, not a resumable campaign.

## Plots

The comprehensive suite adds outputs 70--74: mass-plane status, a 4--9.2 GeV
zoom, maximum bound ratio, two lepton-product panels with individual limit
curves, and mass versus mixing with the invisible threshold marked. Zero
predictions/mixings omitted from logarithmic axes are counted. The ratio colour
scale is linear below one and logarithmic above it, retaining zero signals.

Status maps are generated even for uniform or unassessed samples. Missing
diagnostic plots have explicit omission reasons in the index and summary.
Legacy files have flavour-unassessed status, zero current full-viability
survivors, and a separate labelled `pre_flavour_full_viability` count. The index
reports otherwise viable points removed by flavour and those still unassessed.
The historical collaborator cumulative diagnostic remains its named
HB/HS/W-mass/DM sequence; it is not the full-viability selection.
