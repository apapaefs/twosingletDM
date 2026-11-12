Two-singlet DM scan: odysseus-100x1000-66662
Handover instructions, 29 October 2026

Files to accompany this kit
--------------------------
combined_points_mg5_hh_gw.tsv (or the same file compressed as .tsv.gz)
combined_points_mg5_hh_gw.metadata.json
campaign_state.json

The TSV is the final combined point table, including the original constraint
assessments and the additional MG5, di-Higgs and CalcGW results. Earlier
combined_points*.tsv files represent earlier stages of the same point sample;
they are not independent scans to concatenate. The metadata records the later
calculations, while campaign_state.json records the original scan setup.
applied-constraints.pdf in this kit explains the implemented constraints.
The reference-seed-66662.metadata.json snapshot in this kit gives the detailed
sampling bounds, fixed parameters and software provenance for the first seed.

The original campaign is on Odysseus at:
/home/apapaefs/Projects/TwoSingletDM/twosingletDM/output/odysseus-100x1000-66662

First look: Python 3 only
------------------------
Unzip this kit and put the shared TSV and metadata beside inspect_points.py.
Run from that directory:

  python3 inspect_points.py combined_points_mg5_hh_gw.tsv.gz

For an uncompressed table, remove .gz from the input filename. The reader
streams the input and uses only the Python standard library. It does not run
HiggsTools, micrOMEGAs, BSMPT or MadGraph. It writes:

  quicklook/summary.json                  counts and the exact selection
  quicklook/fully_viable_points.tsv.gz     complete rows for viable points
  quicklook/fully_viable_overview.csv      selected columns for a first look

Open the small overview CSV in a spreadsheet, pandas or Mathematica. The
compressed viable-point table includes every original column, including the
detailed JSON cells. A second invocation requires a new --out directory.

Expected numbers in the audited 29 October snapshot:
  All evaluated parameter points                    591,340
  Points passing evo and thc                        100,000
  Fully viable points                                   571
  Parent points with successful CalcGW predictions    20,763
  Fully viable parent points with successful CalcGW      431

The campaign target was 100 seeds times 1,000 evo-and-thc passes. The complete
file also includes points failing constraints. It is not a sample of 100,000
fully viable points, and the density of sampled points is not a probability
distribution or likelihood.

Selections and identifiers
--------------------------
Each row is one parameter point, identified by the pair (seed, point_index).
point_index alone is not unique across seeds.

For trsm_constraints_v2, the current plotting selection is:

  full viability = thc & experimental_subset & flavour & dm

All four must explicitly be True. Non-DM viability drops only the dm test.
Use experimental_subset for the combined experimental gate; it handles the
applicable coverage of the separate W-mass test. Outside-coverage non-vetoes
are not measurements validating the point. The individual assessment/status
columns give the details.

evo, vacuum_tree_global, rg_bfb, rg_unitarity and theory_strict_subset are
separate diagnostics; they are not extra conditions in this full-viability
definition. In particular, the evo-and-thc stopping rule is a different
selection. See the constraint PDF before adding other theory requirements.

True/False are literal text values. nan or an empty cell means unavailable or
unassessed, not a physical zero. Do not turn strings into Boolean values with
bool(value) or astype(bool): the nonempty string 'False' would be truthy.
For *_excluded fields, True means excluded; this is the opposite sense to the
aggregate pass flags above. Missing cross sections or failed GW jobs are not
zero signals.

Useful columns and conventions
------------------------------
M2, M3                 masses in GeV; M3 is the stable DM scalar mass
vs                     singlet vacuum expectation value, GeV
a12                    visible-scalar mixing angle, radians
k1, k2                 visible coupling modifiers
lX, lPhiX, lSX         dimensionless quartic couplings
K133, K233             potential coefficients of h1*h3^2 and h2*h3^2, GeV
w1, w2                 total scalar widths, GeV
*_br                   dimensionless branching fractions
dm_omega               Omega h^2
dm_dir_det             spin-independent direct-detection cross section, pb
dm_dir_det_limit       comparison limit, pb, already adjusted for abundance
dm_indirect_ratio      implemented gamma-line signal/limit ratio
dm_cmb_ratio           implemented abundance-adjusted CMB signal/limit ratio
mg5_xsec_*_pb          production cross sections in pb, at 13.6 TeV here
mono_higgs_xsec_pb      production rate times the specified invisible BR, pb
mono_z_total_xsec_pb   total ZH2 production rate times BR(H2 -> H3 H3), pb

M1 is fixed at 125.09 GeV. M2 and M3 are independently scanned, with M2 in
[4,1000] GeV and M3 in [4,65] GeV; the labels do not impose ascending masses.
vx=0 and a13=a23=0. K133/K233 are potential coefficients, with corresponding
vertices -2i K133/-2i K233. Other trilinear conventions need the same check.

Simple plot of the viable points
-------------------------------
If pandas and matplotlib are installed:

  import pandas as pd
  import matplotlib.pyplot as plt
  d = pd.read_csv('quicklook/fully_viable_overview.csv')
  p = plt.scatter(d.M2, d.M3, c=d.dm_omega, s=12)
  plt.colorbar(p, label=r'$\Omega h^2$')
  plt.xlabel(r'$M_2$ [GeV]')
  plt.ylabel(r'$M_3$ [GeV]')
  plt.tight_layout()
  plt.savefig('viable_mass_plane.png', dpi=180)

To compare excluded and viable points, load selected columns from the full
file with pandas.read_csv(..., sep='\t', usecols=[...], chunksize=50000).
Pandas reads .gz files directly. Avoid loading all the long JSON columns when
only masses, couplings and constraint flags are needed.

Thermal and GW results
----------------------
ewpt_ew_entry_jump_over_T is the critical EW-entry diagnostic used for the
stars in plot 34b. Stars require a successful stored entry assessment and a
ratio strictly greater than one. It is not a LISA SNR or an extra viability
cut.

gw_transitions is a JSON list. Each item describes one coexistence-pair
transition, so a parameter point can contribute several GW markers. Use
items with valid=true. Fields include transition_index, T_star_GeV, alpha,
beta_over_H and snr_lisa (the native three-year LISA SNR). Join results using
(seed, point_index, transition_index). A star labels the parent parameter
point; it does not assert that every associated GW transition is EW entry.

The all-points CalcGW view does not impose DM viability. For the viable GW
subset, apply the same four pass flags as above before reading transitions.
The SNR=5 line is a reference for interpretation, not part of the stored
viability selection. For a specific interesting point, quote its seed and
point_index so the original detailed output can be located.

Message to send with the files
-----------------------------
Here is the combined table from our low-DM-mass scan, together with its
metadata and the constraint summary. It contains all 591,340 evaluated
parameter points, including excluded points. The current full-viability
selection leaves 571 points. The supplied reader extracts these into a small
CSV for a first look and a complete TSV for detailed inspection. You only
need Python 3 to run it. Please identify interesting points by the pair
(seed, point_index). The later collider and GW results are already included;
missing entries mean unavailable calculations, not zero signals.
