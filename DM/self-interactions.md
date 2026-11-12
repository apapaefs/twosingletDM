# Dark-matter self-interactions (H3 H3 -> H3 H3)

Three separate codes:

- `trsm_sidm.py` computes the H3 self-scattering cross section σ(v) and the
  viscosity cross section σ_V(v) of a vx=0 TRSM point, at fixed relative
  velocities and as Maxwell-Boltzmann halo averages. NumPy only.
- `trsm_sidm_bounds.py` compares σ_V with a table of astrophysical upper
  limits on σ/m (`datafiles/sidm_limits_v1.json`) and writes a verdict.
- `trsm_sidm_lower_limit.py` tests the lower limit σ/m ≳ 1 cm²/g at galactic
  velocities (`datafiles/sidm_lower_limits_v1.json`), the condition for
  self-interactions to produce galactic cores.

`plot_trsm_sidm.py` makes the figures (needs matplotlib). `DM/main_sidm.c`
and `trsm_sidm_calchep_check.py` cross-check the amplitude against CalcHEP.

## Cross section

### Conventions

Canonical-v1 potential coefficients,
`V ⊃ K133 h1 h3² + K233 h2 h3² + (λX/4) h3⁴` with `h3 = X`. The Feynman
rules are −2iK_i33 and −6iλX, as in `models/h4GOn/lgrng1.mdl`. Contact, s-,
t- and u-channel exchange of h1 and h2 contribute:

    M = −6λX − 4 Σ_i K_i33² [1/D_i(s) + 1/(t − M_i²) + 1/(u − M_i²)]

v is the relative velocity (v = 2β_cm, s = 4M3²/(1 − v²/4), k ≈ M3 v/2). For
identical particles

    dσ/dΩ = |M|²/(64π² s),   σ = ½ ∫dΩ dσ/dΩ,   σ_V = (3/2)·½ ∫dΩ sin²θ dσ/dΩ

so σ_V = σ for isotropic scattering. Written with dσ/dcosθ normalised to
σ = ∫dcosθ dσ/dcosθ over cosθ from −1 to 1, this is
σ_V = ∫dcosθ (3/2) sin²θ dσ/dcosθ: the weight of Yang & Yu (arXiv:2205.03392,
Eq. 2.5) and of Ando et al. (arXiv:2403.16633, Eq. 2.25). The σ_V of
Colquhoun et al. (arXiv:2011.04679, Eq. 28) is 2/3 of it.

The exchange terms are kept exactly:
1/(t − M_i²) = −(1/M_i²)/(1 + r_i sin²(θ/2)) with r_i = (M3 v/M_i)²/(1 − v²/4),
and the same with cos²(θ/2) for the u-channel. No expansion in r_i is made.
The t-channel term alone is the angular dependence of Ando et al., Eq. 2.24,
with w = M_i/M3. Halo averages use f(v) ∝ v² exp(−v²/(4σ_1D²)), whose
mean relative velocity is ⟨v⟩ = 4σ_1D/√π, and report ⟨σ_V vⁿ⟩/⟨vⁿ⟩. n = 1 is
the ⟨σv⟩/⟨v⟩ of Sagunski et al. (arXiv:2006.12515). n = 5 is the effective
cross section of Yang & Yu.

### What the default mode does

1. **s-channel propagator.**
   `D_i(s) = s − M_i² + i M_i Γ_i,other + Π_i(s) − Re Π_i(M_i²)`, with the
   H3-loop self-energy
   `Π_i(s) = K_i33²/(8π²) [−β ln((1+β)/(1−β)) + iπβ]`. The imaginary part is
   the energy-dependent h_i → h3h3 width K_i33² β(s)/(8π). The real part is
   subtracted on shell, so M_i is the physical mass. For M_i < 2M3 the
   subtraction leaves `+K_i33² β̄_i/(8π)` at threshold, with
   β̄_i = √(4M3²/M_i² − 1). This is the continuation of the threshold cusp
   below threshold. It keeps the pole of the amplitude at s = M_i².
2. **Unitarity.** Every even partial wave is unitarised,
   `a_l = K_l/(1 − i(β/2)K_l)`, so σ_l ≤ 8π(2l+1)/k².
3. **Non-perturbative exchange.** t/u-channel h_i exchange between slow
   H3's is an attractive Yukawa potential with α_i = K_i33²/(4π M3²). The
   Born approximation needs b_i = α_i M3/M_i ≪ 1. Its relative error on the
   scattering length is b/2, and bound states appear at b = 1.68 (s-wave)
   and 21.9 (d-wave). When max b_i ≥ 0.01 (`--yukawa-threshold`) the Born
   partial waves of the exchange are replaced by phase shifts from the
   radial Schrödinger equation:
   `K_l = K_l^tree − K_l^Born-exchange + (2/β) tan δ_l`.
   Both mediators enter the potential. The contact and s-channel terms are
   kept additively in the K-matrix. Their distortion by the potential is
   not included.

Angular integrals use nodes that cluster towards the forward and backward
directions when k exceeds the lightest mediator mass. They agree with the
closed form to 10⁻⁸ up to k = 1000 M_i.

### Regimes

The column `sidm_regime` says how far the calculation can be trusted.

| Regime | Condition | Treatment |
|---|---|---|
| `perturbative` | max b_i < 0.01 | Born exchange, unitarised |
| `nonperturbative` | max b_i ≥ 0.01 and max α_i < 0.5 | Schrödinger phase shifts |
| `strong` | max α_i ≥ 0.5 (K_i33 ≥ 2.5 M3) | Numbers written, `sidm_reliable = False` |

In the strong regime bound-state velocities are of order α/2, and for
M_i ≳ M3 the exchange is shorter-ranged than the H3 Compton wavelength. A
non-relativistic potential is then not a controlled approximation. The
threshold is adjustable with `--alpha-strong`. A point with vs = 500 GeV,
λSX = 1/2 (K233 = 250 GeV) is in this regime for M3 < 100 GeV.

Other approximations:

- The annihilation channel H3 H3 → h2 h2, open for M3 > M2, is not included
  in the unitarisation. Where it is large the relic abundance is small.
- Loop corrections other than the H3 bubble in the s-channel are neglected.

### What changed in v2 (`sidm_method = trsm_sidm_v2`)

A review of v1 found three problems.

- **Threshold cusp missing.** v1 kept only the imaginary part of the H3
  loop. For M_i slightly below 2M3 this misplaced the pole and overestimated
  σ by (1 + K²/(32π M3 κ))², where κ is the binding momentum. The factor is
  79 for M3 = 10 GeV, M2 = 2M3 − 1 keV, K233 = 5 GeV. Above threshold
  nothing changes.
- **Exchange only flagged.** v1 kept the Born exchange at any coupling and
  set a flag above b = 1. For M3 = 60, M2 = 14, K233 = 120 GeV (b = 1.36)
  v2 is 21 times larger than v1.
- **Angular quadrature.** v1 used fixed Gauss nodes in cos θ and lost the
  forward peak for k ≫ M_i (a factor 5 at k = 100 M_i). The scan ranges
  (M2 ≥ 4 GeV) do not reach that region.

On the 571 viable points of `micromegas-fully-viable-571`, the 519
perturbative points change by less than 1%. Five of the 23 non-perturbative
points change by more than 1% (factors 0.7 to 1.7), and the 29 strong
points change by up to three orders of magnitude.

### Usage

```bash
# one point (widths: --w1/--w2 physical, or --gamma1-other/--gamma2-other)
python trsm_sidm.py point --m2 100.0000125 --m3 50 --vs 500 --a12 0 \
    --lx 0.1 --lphix 0 --lsx 0.004 --gamma2-other 1e-9

# sigma/m and sigma_V/m versus v
python trsm_sidm.py curve ...same point options... --vmin 1 --vmax 1e4 --n 400

# add SIDM columns to a scan (new file plus <output>.sidm.json; .gz accepted)
python trsm_sidm.py scan output/trsm_points_X.dat --output output/trsm_points_X_sidm.dat
```

Scan rows need M2, M3, lX, K133/K233 (or lPhiX, lSX, vs, a12) and the
physical widths w1, w2. The widths are read from the row, not recalculated:
Γ_i,other = w_i − `h_i_h3h3_width`, with the saved `h1_h3h3_width` and
`h2_h3h3_width` (recomputed from K_i33 only if the columns are absent). The
one width the code computes itself is the energy-dependent h_i → h3h3
width K_i33² β(s)/(8π√s), which is needed away from s = M_i² and equals the
saved column on shell. Using the micrOMEGAs widths `dm_h1_width_GeV` and
`dm_h2_width_GeV` in place of w1, w2 changes σ_V by at most 0.2% on the 571
viable points. A perturbative point takes about 0.02 s
and a non-perturbative one up to 1 s. Output columns (σ/m in cm²/g,
velocities in km/s):

- `sidm_sigma_over_m_v<V>kms`, `sidm_sigmaV_over_m_v<V>kms` at fixed v.
- `sidm_sigmaV_over_m_avg_n<n>_disp<D>kms` for 1D dispersion D.
- `sidm_regime`, `sidm_reliable`, `sidm_alpha_h1/h2`, `sidm_born_h1/h2`.
- `sidm_v_res_h1/h2_kms`, `sidm_dv_res_h1/h2_kms`: velocity and width of an
  s-channel pole above threshold.
- `sidm_tree_unitarity_max`: above 1 means the tree amplitude violates
  unitarity at one of the fixed velocities.
- `sidm_avg_quadrature_relerr`, `sidm_status`, `sidm_error`.

## Bounds

```bash
python trsm_sidm_bounds.py list --set all --notes
python trsm_sidm_bounds.py scan output/trsm_points_X.dat --output output/trsm_points_X_sidm_bounds.dat
python trsm_sidm_bounds.py point ...point options... --omega 0.12
```

### Limit table

`datafiles/sidm_limits_v1.json` holds each limit with its confidence level,
the velocity it refers to, the reference and notes. The numbers were
checked against the arXiv abstracts and texts in November 2026. All limits
were derived for velocity-independent, isotropic scattering.

| id | σ/m [cm²/g] | Evaluation | Set |
|---|---|---|---|
| `groups_sagunski2021` | < 1.1 (95%) | ⟨v⟩ = 1150 km/s | default |
| `clusters_sagunski2021` | < 0.35 (95%) | ⟨v⟩ = 1900 km/s | default |
| `mergers_wittman2018` | < 2 (95%) | min over 500–4000 km/s (velocity assumed) | default |
| `clusters_andrade2022` | < 0.13 (95%) | ⟨v⟩ = 1458 km/s | stringent |
| `clusters_eckert2022` | < 0.19 (95%) | ⟨v⟩ = 1000 km/s (velocity assumed) | stringent |
| `bcg_wobble_harvey2019` | < 0.39 (95%) | ⟨v⟩ = 1500 km/s (velocity assumed) | stringent |
| `bullet_randall2008` | < 1.25 (68%) | min over 3100–4700 km/s | stringent |
| `draco_read2018` | < 0.57 (99%) | ⟨v⟩ = 58 km/s (velocity assumed) | dwarf |
| `segue1_hayashi2021` | < 0.086 (2σ) | ⟨v⟩ = 10 km/s (velocity assumed) | dwarf |
| `willman1_hayashi2021` | < 0.39 (2σ) | ⟨v⟩ = 10 km/s (velocity assumed) | dwarf |

`stringent` contains `default`. `all` contains everything. The dwarf limits
are in no other set: they depend on the halo modelling and conflict with
analyses that prefer 20–150 cm²/g for the same galaxies (Correa,
arXiv:2007.02958).

For a Maxwellian entry the model value is ⟨σ_V vⁿ⟩/⟨vⁿ⟩ with
σ_1D = √π⟨v⟩/4 and n = 1 (`--moment 5` gives the Yang & Yu effective cross
section). For a velocity range it is the smallest σ_V in the range.
"Velocity assumed" means the paper states no usable velocity; the notes in
the table say where the number comes from. Robertson et al.
(arXiv:1605.04307) consider the Bullet Cluster limit overstated, and
Wittman et al. advise caution with theirs.

### Sub-dominant H3

With ξ = min(1, Ωh²/0.12) from `dm_omega`, each ratio is
`ξ^p × value / limit` with p = 2 (`--abundance-power`). No bound is applied
for ξ < 0.1 (`--min-fraction`). `--no-rescale` treats H3 as all of the dark
matter.

The limits constrain the total dark-matter distribution. Only a fraction ξ
of the mass can scatter, and each H3 particle sees a target density reduced
by ξ. The number of scatterings per unit total mass is therefore ξ² ρ σ/m,
which gives p = 2 as long as each H3 particle scatters at most about once.
For larger cross sections the effect saturates, because at most a fraction
ξ of the mass is affected. This is the reason for the cut at ξ < 0.1: a
component of 10% or less can have arbitrarily strong self-interactions
(Pollack et al., arXiv:1501.00017; Fan et al., arXiv:1303.1521), and the
Bullet Cluster allows 20–30% (Fan et al.; Choquette et al.,
arXiv:1812.05088). p = 1 rescales only the scattering probability of an H3
particle and overstates the constraint. None of these papers gives a
rescaling law for cluster limits, so p and the cut remain choices.

### Output

- `sidm_<id>_sigma_over_m`: model value before rescaling.
- `sidm_<id>_ratio`: rescaled value over the limit.
- `sidm_bounds_max_ratio`, `sidm_bounds_strongest`.
- `sidm_abundance_fraction`, `sidm_bounds_applicable` (ξ ≥ 0.1).
- `sidm_bounds_unitarity_ratio`: the largest rescaled ratio possible if all
  even partial waves up to l = 8 saturate unitarity. It needs no dynamics.
- `sidm_excluded`: `True`, `False`, or `nan` when unavailable.
- `sidm_bounds_status`, `sidm_bounds_reason`.

`sidm_excluded` is `nan` when the calculation failed, the abundance is
missing, or the point is in the strong regime while the bounds would apply.
A strong point still gets `False` if `sidm_bounds_unitarity_ratio < 1` and
the typical momentum is below the lightest mediator mass. Unavailable must
not be read as passing.

### Result for the 571 viable points

| Run | Excluded | Not excluded | Unavailable |
|---|---|---|---|
| default set | 0 | 571 | 0 |
| all limits | 0 | 568 | 3 |
| all limits, `--no-rescale` | 0 | 542 | 29 |

400 points have ξ < 0.1. The largest model value is 3×10⁻³ cm²/g and the
largest rescaled ratio 1×10⁻⁷. The unavailable points are strong-coupling
points for which the dwarf limits are not out of reach by unitarity.

Unitarity explains the gap. One partial wave gives at most
σ/m = 8π(2l+1)/(k² M3), which for the s-wave is about 2000/M3³ cm²/g at
1000 km/s (M3 in GeV) and 2×10⁶/M3³ cm²/g at 30 km/s. An s-wave resonance
can therefore reach the default limits only for M3 below 13 GeV, the
stringent set below 24 GeV, and the Draco limit below 110 GeV. The 571
points have M3 between 12 and 65 GeV.

## Lower limit (galactic cores)

```bash
python trsm_sidm_lower_limit.py list --set all --notes
python trsm_sidm_lower_limit.py scan output/trsm_points_X.dat --output output/trsm_points_X_sidm_lower.dat
python trsm_sidm_lower_limit.py point ...point options... --omega 0.12
```

This is a requirement, not an exclusion. A point that fails it is allowed;
its dark matter is collisionless in galaxies and does not produce cores.
It contradicts the `dwarf` set of upper limits (for example Draco, below
0.57 cm²/g at 58 km/s), so the two cannot both be imposed.

### Table

`datafiles/sidm_lower_limits_v1.json`, checked against the papers in
November 2026:

| id | σ/m [cm²/g] | Mean relative velocity | Set |
|---|---|---|---|
| `cores_tulin_yu2018` | ≥ 1 | 30–200 km/s | default |
| `rotation_curves_ren2019` | ≥ 1 | 29–430 km/s (velocity assumed) | all |
| `tbtf_milky_way_tulin_yu2018` | ≥ 0.6 | 50 km/s | all |
| `tbtf_local_group_tulin_yu2018` | ≥ 0.5 | 50 km/s | all |

- Tulin & Yu (arXiv:1705.02358, Table I) give σ/m ≳ 1 cm²/g at
  v_rel = 30–200 km/s for cores in dwarf and LSB spirals, and the
  too-big-to-fail values.
- Ren et al. (arXiv:1808.05695) fix σ/m = 3 cm²/g for 135 SPARC galaxies
  with V_f from 20 to 300 km/s and state that results are the same above
  about 1 cm²/g. Kamada et al. (arXiv:1611.02716) fix 3 cm²/g for 30
  galaxies with V_f from 25 to 300 km/s. Neither fits a lower limit.
- Ren et al. use ⟨σ v_rel⟩ = σ (4/√π) σ_v0, with σ_v0 the 1D dark-matter
  dispersion. The range 29–430 km/s assumes σ_v0 = 0.64 V_max (Yang & Yu)
  and V_max = V_f, so ⟨v⟩ = 1.44 V_f.

### Evaluation

The model value is ⟨σ_V v⟩/⟨v⟩ on a logarithmic grid of mean relative
velocities across the range. The smallest value is compared with the
minimum, because the requirement has to hold in all those halos.

The value is multiplied by ξ² as for the upper limits. The requirement is
not met for ξ < 0.9 (`--min-fraction`), whatever the cross section: the
collisionless remainder of the dark matter keeps its cusp.

Output: `sidm_lower_<id>_min`, `_max` (before rescaling) and `_ratio`,
`sidm_lower_min_ratio`, `sidm_lower_weakest`,
`sidm_lower_unitarity_ratio`, and the verdict `sidm_lower_satisfied`
(`True`, `False`, or `nan` when unavailable).

### Result for the 571 viable points

None satisfies the lower limit.

- 13 points have ξ ≥ 0.9. Their largest σ_V/m is 1×10⁻⁸ cm²/g.
- The largest value among reliable points is 4×10⁻⁶ cm²/g, and the largest
  of all is 3×10⁻³ cm²/g (a strong-coupling point with ξ = 3×10⁻¹⁰).
- No point has a velocity-dependent cross section: the values at 30 and
  200 km/s differ by less than 0.3%.

What it would take: the s-wave unitarity ceiling at ⟨v⟩ = 200 km/s exceeds
1 cm²/g only for M3 below 40 GeV, and reaching it needs a bound state at
threshold. An example is M3 = 20 GeV, K233 = 5 GeV and
2M3 − M2 = 0.8 keV. It gives 9 cm²/g at 30 km/s and 2 cm²/g at 200 km/s,
and stays below all cluster limits because it falls like 1/v² above. The
mass tuning is 2 parts in 10⁸.

### Can a thermal H3 satisfy it?

The cross section can be reached, but not together with the relic
abundance and the CMB bound. Both routes were run through micrOMEGAs 7.1.4
with K133 = 0.

- **h2 pole at threshold** (the example above). sin θ = 1.5×10⁻³ gives
  Ωh² = 0.121 and the lower limit is satisfied. The same pole enhances
  s-wave annihilation into SM particles at low velocity: σv today is
  3×10⁻¹⁹ cm³/s and the Planck CMB ratio is 1×10⁷. For sin θ from 10⁻⁵ to
  1.5×10⁻³ the abundance is at or above 0.12 and the CMB ratio is 450 or
  more. For larger sin θ the abundance drops below 0.9 of the total.
- **Bound state of the exchange potential at threshold** (M3 > M2,
  b = 1.68). M3 = 20, M2 = 5, K233 = 45.9 GeV gives 7.9 cm²/g, but
  Ωh² = 8×10⁻⁸ because H3 annihilates into h2 pairs. The coupling that
  gives Ωh² = 0.12 there is K233 = 1.3 GeV, for which σ_V/m is
  6×10⁻¹⁴ cm²/g.

micrOMEGAs uses fixed widths and no threshold term, so the numbers are
indicative; the gaps are many orders of magnitude. An abundance of
non-thermal origin with sin θ below about 5×10⁻⁷ would pass the CMB bound
in the first case. The scan assumes thermal freeze-out and does not cover
that.

## Figures

```bash
python plot_trsm_sidm.py landscape --m3 60 --k233 1 10 50 120 --output-dir plots/sidm
python plot_trsm_sidm.py velocity --output-dir plots/sidm
python plot_trsm_sidm.py scan scan.dat --output-dir plots/sidm
python plot_trsm_sidm.py points scan.dat --output-dir plots/sidm
```

- `landscape`: σ_V/m against M2/M3 from 0.07 to 3.5, covering M3 > M2 and
  the s-channel pole at M2 = 2M3.
- `velocity`: σ_V/m against v for M3 > M2 benchmarks, with the unitarised
  Born result for comparison.
- `scan`: the points of a scan against two limits, marked by regime and by
  whether ξ ≥ 0.1.
- `points`: ⟨σ_V v⟩/⟨v⟩ against ⟨v⟩ for a few points of a scan (the largest
  cross sections, the points closest to the two s-channel poles, the
  largest with M3 > M2), with the lower limit, the group and cluster upper
  limits and the tuned example above.

For M2 of a few GeV or more, σ_V is independent of velocity up to 10⁴ km/s
unless a bound state sits at threshold. For M3 = 60 GeV and K233 up to
120 GeV, σ_V/m exceeds 10⁻⁴ cm²/g only close to the s-channel pole or to a
Yukawa bound-state threshold.

## Validation

```bash
python -m unittest test_trsm_sidm.py test_trsm_sidm_bounds.py test_trsm_sidm_lower_limit.py

cp DM/main_sidm.c ../runtime-v2/micromegas_7.1.4/TRSM/
make -C ../runtime-v2/micromegas_7.1.4/TRSM main=main_sidm.c
python trsm_sidm_calchep_check.py --main-sidm ../runtime-v2/micromegas_7.1.4/TRSM/main_sidm
```

- **CalcHEP.** `main_sidm.c` evaluates `~X,~X->~X,~X` with the model's own
  vertices and widths. CalcHEP keeps a propagator width only within
  |q| ≲ 2.7 MΓ and, with a large range, also in t-channel propagators
  (`sqme_aux.c:prepDen`). With that prescription the trsm_sidm tree
  amplitude reproduces CalcHEP totals and cos θ bins to ≤ 8×10⁻⁷ for seven
  benchmarks at v/c from 10⁻⁴ to 0.5 (micrOMEGAs 7.1.4). CalcHEP uses fixed
  widths and no unitarisation, so it differs from the default mode near
  M2 ≈ 2M3.
- **Self-energy.** Tests check the dispersion relation between the real and
  imaginary parts, the pole position for M_i < 2M3, and the limit
  σ → 8π/(κ² + k²) of a strongly coupled shallow pole.
- **Schrödinger solver.** Phase shifts agree with an independent SciPy
  integration to 10⁻⁶. Halving the step changes the phase by about 2×10⁻⁷
  at most, for b from 0.01 to 250, l ≤ 8 and k up to 20 M_i. The bound-state thresholds
  come out at b = 1.6798 and 21.895.
- **Viscosity.** σ_V reproduces 3/2 times the phase-shift formula of
  Colquhoun et al. for identical bosons. Integrating Eq. 2.24 of Ando et
  al. with the weight (3/2) sin²θ reproduces their Eq. 2.25.
- **Exchange at large momentum.** The tree cross sections agree with an
  independent integration to 10⁻⁹ for r_i from 10⁻⁴ to 10⁴.
