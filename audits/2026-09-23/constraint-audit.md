# Constraint audit — 23 September 2026

## Scope and conclusion

Audited repository commit **`b7b0d934fe85a45842174af28cf03005f41fbda4`**, present on local `main` and on manto. This audit covers the active `vx=0` scan, its HiggsTools, theory, EW precision and DM gates, BSMPT integration, phase-transition labels, and the new freeze-out diagnostics. It includes source inspection, primary documentation, numerical counterexamples, and regression tests in manto's installed environment.

**The agreed baryogenesis and GW candidate formulas are implemented correctly. Several other flags need qualification or repair before “passes all constraints” can be interpreted as a robust physical exclusion/acceptance statement.** The most consequential findings are inconsistent SM inputs between programs, absence of a global-vacuum requirement in the main selection, the weak meaning of the RG flag, and the effective photon/gluon vertices in the DM model. BSMPT status handling and the reconstructed thermal history also need attention.

This is an audit of the current implementation, not a reevaluation of every historical scan point. Only audit artifacts were added; constraint code, running campaigns and external installations were not modified. Priorities below identify scientific impact and repair order; they do not claim that every stored point is affected.

## 1. What the scan actually selects

For the active dark-matter branch,

\[
P_{\rm nonDM}=\texttt{evo}\land\texttt{thc}\land\texttt{hb}\land
\texttt{hs}\land\texttt{ewpo}\land\texttt{wmass},\qquad
P_{\rm full}=P_{\rm nonDM}\land\texttt{dm}.
\]

Neither EWPT candidate flag, the freeze-out Z2 flag, nor the thermal-input sensitivity flag enters `P_full`. A point may be saved despite failing it, depending on output options. The implementation is in [generate_trsm_points.py](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/generate_trsm_points.py:1635).

| Item | Current numerical rule | Effect |
|---|---|---|
| HiggsBounds `hb` | `bounds(pred).allowed` | Required by `P_full` |
| HiggsSignals `hs` | `chi2(point) - chi2(SM) < 4.00` | Required |
| Tree theory `thc` | Quartic copositivity and high-energy scalar unitarity | Required |
| Evolution `evo` | Successful one-loop ODE integration from 91 to 1000 GeV; finite saved solution at 900 GeV | Required |
| EW precision `ewpo` | Correlated STU chi-square `<= 7.82` | Required |
| W mass `wmass` | Tabulated bound on `abs(sin(a12))`, for `133 <= M2 <= 999` GeV | Required, but automatically passes outside coverage |
| Relic abundance | `Omega h² <= 0.121`; no lower bound | Part of `dm` |
| Direct detection | `xi * sigma_SI(neutron) <= limit(M3)`, `xi=min(1,Omega/0.121)` | Part of `dm` |
| Fermi gamma lines | Each available line flux, rescaled by `xi²`, below its R16 bound | Part of `dm`; unavailable lines do not exclude |
| CMB, when enabled | `min(1,Omega/0.12)² * R_Planck <= 1` | Part of `dm`; unavailable enabled CMB result fails aggregate DM |
| Baryogenesis candidate | EW-symmetric to EW-broken FOPT with critical EW jump/T strictly `>1` | Annotation |
| GW candidate | Any-field jump/T strictly `>1` at any available `Tc`, `Tn`, or `Tp` | Annotation |
| Z2 compatibility | No sampled global X-breaking at or below nominal freeze-out temperature, with required coverage | Annotation |
| Thermal sensitivity | At least 10% EW or S VEV change relative to T=0 in `[Tf/2,2Tf]` | Annotation |

### Execution and retention matter

- `--run-ewpt` runs BSMPT only for fully accepted points.
- To retain the collaborator's interest in DM-failing points, use **both** `--run-ewpt` and `--run-ewpt-on-dm-failed`. The latter covers points passing all non-DM gates and creates the DM-failed output. `--write-dm-failed-to-main` additionally retains those rows in the main output.
- `--write-evo-thc-points` saves theory/evolution-passing rows even if experimental or DM flags fail. It does **not** by itself run BSMPT on all of those rows.
- `--write-all-points` runs BSMPT for every evaluated point and bypasses the optional Eq. 4.18 prefilter.
- EWPT reprocessing selects `evo & thc & hb & hs & ewpo & wmass`, independently of DM.
- The current command-line entry point uses the `vx=0` branch. The old general three-VEV evaluation function remains in the source but is incomplete: its theory/RG calls are commented out, it sets both flags to `True`, and the HiggsTools adapter does not populate H3. It must not be reused as a validated general-TRSM constraint path.

### Sampling is not an exclusion constraint

The default mass sampler imposes a conditional lower endpoint `M3 >= M2 + mhiggs`; this is a sampling choice. `--approximate-resonantDM` still samples both `M2 ≈ 2 M3` and `M3 ≈ 2 M2`. In the Z2-odd-X model, the s-channel annihilation poles are `2 M3 ≈ M1` and `2 M3 ≈ M2`. The reverse relation `M3 ≈ 2 M2` is not such a pole. Exact `--resonantDM1` uses `M3=125.09/2`, which currently disagrees with the mass used by micrOMEGAs (finding A).

## 2. Findings requiring attention

### A — High priority: the programs do not use the same SM point

The micrOMEGAs card writes only seven BSM inputs; it does not override `Mh`, `EE`, `MW`, or `SW`. Both installed versions on manto retain the model-file defaults:

| Input/derived value | Scan/HiggsTools | BSMPT TRSM | micrOMEGAs TRSM |
|---|---:|---:|---:|
| SM-like Higgs mass [GeV] | 125.09 | 125.09 | **125.0** |
| EW VEV [GeV] | 246.0 | 246.2196508 | **243.2099703** |

The micrOMEGAs value is `2*MW*SW/EE`, with `MW=80.385`, `SW=0.474`, `EE=0.31333`. Its derived Fermi constant is `1.1954248e-5 GeV^-2`, compared with BSMPT's `1.1663787e-5` input.

**Impact:** the h1 pole is displaced by 90 MeV in mediator mass, or 45 MeV in the DM mass at half the pole. Thus the scan's “exact h1 resonance” and the recorded h1 resonance gap are not the resonance seen by micrOMEGAs. The EW-VEV mismatch also changes portal trilinears, including cancellations in K133/K233. The smaller 246 versus 246.22 mismatch should be aligned as part of the same convention fix.

**Action:** explicitly supply a common SM input convention to every program and record the values actually used. Reevaluate narrow-resonance and cancellation points after alignment; a metadata label alone cannot correct their rates.

Sources: [DM card](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/test_trsm_DM.py:219), [model defaults](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/DM/models/h4GOn/vars1.mdl:4), [model relations](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/DM/models/h4GOn/func1.mdl:25), [resonance diagnostic](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/dm_thermal_relic_diagnostic.py:48). Installed model values were read directly on manto.

### B — High priority: `thc=True` does not establish the desired vacuum

The tree-level test checks quartics, not the depths of stationary points. Positive physical scalar masses establish local stability at the input vacuum; they do not exclude a deeper X-breaking minimum on this partially broken branch.

A reproduced counterexample has

`vs=300, M2=200, M3=50 GeV, a12=0.1, lX=0.1, lPhiX=0.1, lSX=1`.

It gives **`thc=True`, `evo=True`**, while:

| Stationary point `(h,s,x)` [GeV] | Potential [GeV⁴] | Radial Hessian eigenvalues [GeV²] |
|---|---:|---|
| `(246,300,0)` | `-6.1210e8` | `2500,15647.5,40000` |
| `(0,0,674.728)` | **`-5.1815e9`** | `13342.8,91051.6,206758.5` |

Both are local minima; the second is much deeper. This is a counterexample to the meaning of the theory flag, not a claim that this particular example passes the other experimental/DM gates.

BSMPT's default NLO check can reject the intended vacuum within CalcTemps. However, BSMPT is optional and its vacuum status does not enter `P_full`. Hence it does not close this gap in the main selection.

**Action:** add a separate zero-temperature vacuum result. For this polynomial model a tree-level comparison of all physical stationary branches is inexpensive. Keep the loop-level stability result distinct. If metastable input vacua are to be retained, their acceptability needs a lifetime treatment or an explicit “unassessed metastability” label.

Source: [theory test](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/test_trsm_theory_constraints.py:149). The fully broken TRSM's tree-level global-vacuum statement must not be imported into the `vx=0` branch without checking its assumptions; the example above resolves this directly. [TRSM theory reference](https://arxiv.org/html/1908.08554v2#S3.SS1).

### C — High priority: the RG flag is numerical solvability, not running perturbativity/stability

No bound is applied to running quartics or scattering eigenvalues. Boundedness is not rechecked along the trajectory. The only output-scale test is finiteness at 900 GeV, alongside solver success through 1000 GeV.

For `vs=300, M2=200, M3=50 GeV, a12=0, lX=3, lPhiX=lSX=0`, both input theory and `evo` pass. A tighter integration gives

`lambdaX(900 GeV)=13.8639`, `lambdaX(1000 GeV)=16.6334`.

The code's own tree-level perturbative-unitarity test fails at both scales. The evolution flag therefore cannot be described as stability or perturbativity up to 1 TeV, let alone to a higher scale.

There is also a separate initialization error: `mmS2init` omits `M3²`. With the canonical potential it should be `M3² - lPhiX*v²/2 - lSX*vs²/2`. Dimensionful parameters do not feed the current dimensionless beta functions, so this does not directly change the quartic trajectory, although it makes the mass running incorrect and can influence numerical solver behavior.

**Action:** define the intended running-coupling criterion and check it along the trajectory; retain integration failure as a distinct numerical status. Correct the X mass initialization. Specify the starting matching scale and the treatment of thresholds if a UV-validity claim is intended.

Sources: [integration acceptance](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/test_trsm_evolution.py:22), [mass initialization and scales](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/test_trsm_evolution.py:222).

### D — High priority: DM photon/gluon effective vertices have an extra mixing factor

The model defines `LAAh=-abs(lAAhiggs(Mh,"h1"))` and similarly `LGGh`. The installed CalcHEP implementation obtains these functions from the actual h1 couplings to charged/colored particles. Those h1 couplings already contain `CosT`. The model then multiplies the effective h1 vertex by `CosT` again, and the effective h2 vertex by `SinT`.

At a fixed loop scale, the resulting mixing dependence is proportional to `cos²(a12)` for h1 and `sin(a12)*cos(a12)` for h2, instead of one respective doublet factor. Both states also share a form factor evaluated at `Mh=125 GeV`; its complex phase is discarded with `cabs`. A general annihilation amplitude at invariant mass `sqrt(s)` is not represented by that fixed form factor.

**Evidence level:** the extra factor follows from the installed source and model vertices. This audit has not recomputed a corrected grid of annihilation rates, so its impact on total relic abundance is not quantified. It directly concerns gamma-line predictions and the gluon channel, and may affect widths used in resonant annihilation.

**Action:** validate the effective vertices against the SM and mixed-scalar limits, with one mixing factor and the appropriate loop kinematics, before interpreting those rates quantitatively.

Sources: [loop definitions](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/DM/models/h4GOn/func1.mdl:31), [photon vertices](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/DM/models/h4GOn/lgrng1.mdl:4), [gluon vertices](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/DM/models/h4GOn/lgrng1.mdl:37). Upstream implementation inspected on manto: `/Users/apapaefs/Projects/micromegas_7.1.4/CalcHEP_src/c_source/dynamicME/dynamic_cs.c`, `HiggsLambdas` at line 689 and `lAAhiggs` at 771.

### E — Medium priority: the default direct-detection fit is discontinuous and insufficiently sourced

The inherited `lz2025-source` prescription changes at 10 GeV from approximately `1.844694e-9 pb` immediately below the boundary to `1.094375e-10 pb` at the boundary: a **factor 16.86 tighter** exclusion for an arbitrarily small mass change. Other joins are much smaller. Its low-mass behavior is inherited rather than a validated current LZ mass curve.

The optional normalized-table path is better specified: explicit source, units, observed 90% upper limit, checksum, log-log interpolation, and no extrapolation. The separate LZ 2026 high-mass table is explicitly a temporary approximation, not an official mass-dependent numerical release, and is not the default.

**Action:** replace the default piecewise fit with a verified, documented table or explicitly label the existing choice as an approximation. The experimental comparison cannot be certified from its name alone.

Sources: [fit](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/test_trsm_DM.py:477), [table documentation](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/DM/direct-detection.md:1).

### F — Medium priority: unsuccessful BSMPT physics can become “not a candidate”

`ewpt_status="success"` means the wrapper obtained output without an exception. It does not require `status_nlo_stability`, `status_tracing`, or `status_coex_pairs` to demonstrate a successful physical analysis. With an empty transition list, both new candidate flags become `False`, even when tracing failed or NLO stability prevented the transition calculation.

**Action:** keep the executable result, vacuum status, tracing status and candidate result separate. Use “unknown/unassessed” when a required calculation failed; reserve a negative candidate flag for a sufficiently assessed point. This does not require making nucleation a selection gate.

Sources: [wrapper](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/generate_trsm_points.py:1303), [candidate defaults](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/ewpt_entry_criterion.py:89). The [BSMPT status-code documentation](https://arxiv.org/html/2404.19037v2#S3.SS12) distinguishes NLO stability, tracing coverage and transition-temperature outcomes.

### G — Medium priority: the reconstructed global branch has visible switching artifacts

The Python postprocessor interpolates each phase's VEVs and potential linearly onto the union of its temperature grids, then picks the lowest interpolated potential. It does not reevaluate or minimize the potential at every union-grid temperature.

For archived point 223615 the chosen branch alternates between phase 3 and phase 2 at `50.899,51.000,51.899,52.000,52.899,53.000 GeV`. It reports X-breaking intervals `[51,51]`, `[52,52]`, and `[53,130.779]`. The final-restoration bracket `[50.899,51.000] GeV` does **not** contain CalcTemps' critical crossing `51.5525 GeV`.

This demonstrates an inconsistency between the reconstructed envelope and the resolved transition calculation. Interpolation on offset phase grids is a plausible cause; a direct potential check would be needed to attribute it uniquely. The interval edges must not be treated as certified transition-temperature brackets, particularly for freeze-out near a crossing.

**Action:** reconcile nearby phase switches with the coexisting-phase data, evaluate competing potentials at common temperatures, and retain an unresolved status when the ordering is numerically uncertain. This is distinct from the physical uncertainty caused by supercooling or failed nucleation.

Sources: [interpolation](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/test_trsm_ewpt.py:467), [envelope selection](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/test_trsm_ewpt.py:532), [archived evidence](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/audits/2026-09-23/archived-point-223615.json).

### H — Medium priority: two definitions are still presented as “strong FOPT”

The older plotting classification uses `v_EW,true/T >= 1` for a selected transition, with temperature priority `Tn`, then `Tp`, then completion, then `Tc`. It does not require EW symmetry breaking across that transition. The new baryogenesis classification uses the critical EW-entry jump and strict `>1`.

Point 223615 illustrates the difference: the legacy selection gives `v_EW,true/Tn=1.19153`, but its EW jump is only `5.53e-5`. Its actual EW-entry transition has critical jump `0.13598`, and it correctly fails the new baryogenesis flag.

**Action:** label the old quantity as a legacy true-phase VEV diagnostic and use `ewpt_baryo_candidate` for baryogenesis candidate counts. Preserve the separate GW classification. Mixing these definitions changes the physical interpretation of plots.

Sources: [legacy selection](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/generate_trsm_points.py:1252), [plot classification](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/plot_trsm_constraint_suite.py:1876).

### I — Medium priority: degenerate unitarity eigenvalues can be rejected numerically

The polynomial coefficients match the expected high-energy scalar-scattering cubic. However, `numpy.roots` can produce tiny imaginary components for repeated real roots; the code interprets them as failure.

With `(lambdaPhi,lambdaS,lambdaX)=(125.09²/(2*246²),1,1)` and all portals zero, the exact roots are `(1.55141,6,6)`, comfortably allowed. The implementation produces `6 ± 1.13e-7 i` and returns `False`.

**Action:** use the real symmetric matrix directly,

\[
\begin{pmatrix}
12\lambda_\Phi&2\lambda_{\Phi S}&2\lambda_{\Phi X}\\
2\lambda_{\Phi S}&6\lambda_S&\lambda_{SX}\\
2\lambda_{\Phi X}&\lambda_{SX}&6\lambda_X
\end{pmatrix},
\]

and a symmetric eigensolver. This gives the same intended cubic without the repeated-root instability. The copositivity implementation has the correct half-portal matrix normalization. Its non-strict tolerance also permits quartic-flat boundaries; negative quadratic terms along an exactly flat direction require separate treatment.

Source: [unitarity solver](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/test_trsm_theory_constraints.py:67). The intended polynomial and bounds agree with [TRSM equations 29–33](https://arxiv.org/html/1908.08554v2#S3.SS1).

### J — Smaller numerical and reporting defects

- **HiggsSignals reference offset:** the SM baseline uses HiggsPredictions' built-in SM widths/BRs, whereas the scan overrides them with its external interpolated/normalized arrays. At the exact decoupling limit, `chi2=152.688434` versus baseline `152.548868`, giving a spurious `Delta chi2=0.139565`. Align the reference treatment. This is small relative to the cut of 4, but material for boundary points. [Reference initialization](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/test_trsm_higgstools.py:41).
- **EWPO removable singularity:** `M2=MZ=91.1876 GeV` raises `ZeroDivisionError` in the loop expressions. Nearby values are finite. Implement analytic limits and reject nonfinite computed chi-squares explicitly. [Loop functions](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/singlet_EWPO.py:74).
- **Out-of-range W-mass flags:** outside 133–999 GeV the function warns and returns `True`. Thus a light-scalar point's `wmass=True` does not mean this constraint was evaluated. Save coverage separately. [W-mass test](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/singlet_EWPO.py:45).
- **Relic calculation error and precision:** `darkOmega` supplies `err`, but the one-component driver neither checks nor prints it. `Omega` and `Xf` are printed with only three significant digits; the wrapper makes its cuts on those rounded numbers. For example an underlying abundance 0.1214 would print as 0.121 and pass the upper cut. Preserve full precision and an explicit solver status. [Driver](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/DM/main.c:277).

## 3. HiggsTools and other theory/precision constraints

### What is correctly connected

The active model has `k1=cos(a12)`, `k2=sin(a12)`, `k3=0`. HiggsTools receives masses and scaled SM effective couplings for h1 and h2, physical total widths, visible branching fractions, direct invisible widths, h1→h2h2 and h2→h1h1 when open. The adapter clears point-dependent exotic branching fractions between evaluations. Direct invisible and cascade contributions are included in the common total width, so visible BRs are consistently diluted.

The HiggsBounds verdict is the package's expected-sensitivity-selected exclusion test, using the observed limit of the selected analysis for each particle; it is not a hand-cut on every observed ratio. [Official HiggsBounds API](https://higgsbounds.gitlab.io/higgstools/HiggsBoundsAPI.html).

HiggsSignals supplies the total chi-square. The repository chooses the SM-referenced cut of 4; it is not a profiled confidence region relative to the model's best fit and is not automatically a universal 95% confidence statement for a multi-parameter scan. The saved ranked individual chi-square contributions are useful diagnostics, but their sum need not reproduce the correlated total. [Official HiggsSignals API](https://higgsbounds.gitlab.io/higgstools/HiggsSignalsAPI.html).

### Scope limits

- A provider-only width floor of `1e-10 GeV` avoids HiggsTools discarding very narrow states; the stored physical widths remain unchanged. It does not validate prompt-decay searches for long-lived h2 states.
- The adapter does not supply non-resonant di-Higgs production predictions. HiggsTools emits warnings that these cross sections are not computed. Optional MG5 rates are observables, not additional experimental cuts.
- There is no separate width-to-mass acceptance cut or complete treatment of resonance interference. Applicability of narrow-width experimental limits still needs checking for broad benchmarks. [HiggsTools paper](https://www.sciencedirect.com/science/article/pii/S0010465523001480).
- The tree-level theory criterion contains BFB and asymptotic scalar unitarity; it contains neither finite-energy unitarity nor loop-convergence estimates.
- STU uses `sin²(a12)[O(M2)-O(125.09)]` and a correlated Gaussian fit. The implemented fit is attributed to Jens Erler, private communication, 16 May 2025: central values `(-0.045235,0.013373,-0.016131)`, errors `(0.094703,0.11518,0.087103)`, correlations `(rhoST,rhoSU,rhoTU)=(0.048,0.747,-0.162)`. I could audit the arithmetic, but not authenticate private fit inputs against a public likelihood. The separate W-mass curve and STU fit should not be described as independent pieces of a combined statistical likelihood.

## 4. Dark matter: interpretation of the accepted flag

The wrapper correctly extracts the **neutron SI cross section in pb**, not a signed amplitude. Its DD abundance factor is linear in `xi`; line and CMB signals use `xi²`. Equality with the implemented upper bounds passes. Underabundant points are permitted, implicitly treating X as a component of DM with the assumed corresponding local/cosmological abundance fraction.

The relic upper limit 0.121 is an explicit scan convention, not a relic-density likelihood with theory errors. The CMB normalization instead uses 0.12, as prescribed by that implementation; this small mismatch is transparent but worth keeping explicit.

The default backend is **6.1.15, with CMB disabled**. New 7.1.4 runs enable CMB by default. Legacy resumes retain their saved settings. Enabling CMB requires a compatible driver, checked through `--capabilities`; manto's version-7 driver matches the current `DM/main.c`, while the active version-6 driver predates the CMB additions.

Important physical limitations of the current DM assessment:

1. `darkOmega` is being used in the standard fixed-vacuum freeze-out setup. It does not consume BSMPT's thermal VEVs, phase transition entropy release, or an X-broken cosmological history. Thermal/chemical equilibrium and the assumed expansion history are not independently validated point by point.
2. `VZdecay=VWdecay=0` leaves virtual-vector final-state options disabled. This matters near vector thresholds and for subthreshold mediator decays. [Official micrOMEGAs version documentation](https://micromegasdm.github.io/v6.0/intro_6.0.html).
3. Fermi uses a 2015 R16 monochromatic-line table over 0.214–462 GeV and an approximate ROI integral in 0.5-degree patches. It is not a dwarf-spheroidal continuum exclusion or a general indirect-detection likelihood. Missing/out-of-range/invalid lines are marked unavailable and do not fail DM; hence “passes indirect detection” can mean “no usable line test.” The table's physical context is [the Fermi-LAT line search](https://arxiv.org/abs/1506.00013).
4. CMB uses `PlanckCMB` with the low-velocity spectrum and the s-wave interpretation. The recorded production setting has `vRot=220 km/s`, with the default relative-velocity prescription. The authors explicitly limit this interpretation when the cross section varies at small velocities; it is not a recombination-era velocity-evolution calculation. [micrOMEGAs 7, §6.3.1](https://arxiv.org/html/2606.06645v1#S6.SS3.SSS1).
5. No dedicated DM self-interaction, late-decay/BBN, or collider missing-energy exclusion likelihood is applied by this wrapper. Their absence is a coverage statement, not a proposal to automatically add every possible bound.

### Z2-at-freeze-out and thermal sensitivity

The two diagnostics are appropriately separate. The Z2 flag looks for X-breaking at **any sampled temperature at or below** nominal `Tf`; this also catches later breaking after nominal freeze-out. Missing range/boundary information can remain unknown. The threshold for identifying X-breaking is 1 GeV by default, so “unbroken” is a numerical classification rather than an exact algebraic test of Z2.

The thermal-input diagnostic compares `abs(w1)` and `abs(ws)` against their T=0 values within `[Tf/2,2Tf]`, with a 10% screening threshold. It records boundary overlap separately and avoids VEV interpolation across different identified phases. Resonance fields record signed zero-temperature gaps `Mi-2M3` and their magnitudes divided by `Tf`.

These are sensible **sensitivity screens**, not certificates of relic-density validity. Even small individual VEV changes can strongly affect a cancellation in

\[
K_{133}=\tfrac12(\lambda_{\Phi X}v\cos\alpha-\lambda_{SX}v_s\sin\alpha),\quad
K_{233}=\tfrac12(\lambda_{\Phi X}v\sin\alpha+\lambda_{SX}v_s\cos\alpha).
\]

The screen does not recalculate thermal mixing, propagator masses, widths, couplings or annihilation rates. Gap/Tf does not include the mediator width. If the actual universe remained in a metastable X-breaking phase, an equilibrium-branch `z2_freezeout_compatible=True` cannot establish stability at freeze-out. A restoration only just above freeze-out also does not demonstrate re-equilibration. These limitations are especially relevant to the resonance/cancellation regions motivating this scan.

Sources: [X history](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/ewpt_x_history.py:33), [thermal diagnostic](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/dm_thermal_relic_diagnostic.py:115).

## 5. The precise strong-transition criteria now implemented

### Baryogenesis candidate — matches the agreed critical-only rule

For a CalcTemps critical transition, with the default symmetry tolerance,

\[
|w_{1,f}(T_c)|<5\ \mathrm{GeV},\qquad
|w_{1,t}(T_c)|\geq5\ \mathrm{GeV},\qquad
\frac{|w_{1,t}(T_c)-w_{1,f}(T_c)|}{T_c}>1.
\]

If several critical transitions meet the EW-entry condition, the one with the largest EW jump/Tc is stored. The same transition's Tn and Tp jumps are saved when available. **A missing nucleation, percolation or completion result does not reject this candidate.** A larger nucleation jump cannot rescue a subthreshold critical jump under the agreed rule. Singlet-only and already-EW-broken to EW-broken transitions do not satisfy this EW-entry condition.

This is a qualitative baryogenesis-motivated transition label. Actual preservation of an asymmetry depends on the broken-phase sphaleron washout rate over the relevant history; a fixed one-loop VEV/T threshold is a proxy with gauge/scheme dependence. The present real-scalar extension also does not calculate an additional CP-violating source or a baryon asymmetry. [Baryon-washout analysis](https://arxiv.org/abs/1101.4665).

### GW candidate — matches the agreed any-field rule

\[
\Delta\phi(T)=\sqrt{(\Delta w_1)^2+(\Delta w_x)^2+(\Delta w_s)^2},\qquad
\max_{\substack{\text{reported transitions}\\T\in\{T_c,T_n,T_p\}}}
\frac{\Delta\phi(T)}{T}>1.
\]

Transitions beginning at zero VEV and singlet-only transitions are eligible. A successful Tn or Tp value can qualify a point even if its Tc ratio does not. Approximate nucleation and completion temperatures are saved in the detailed BSMPT output but do not enter this cut. Maxima at different temperature kinds may refer to different transitions. Both criteria use strict `>1`, not `>=1`.

This is a large-jump candidate definition, not a prediction of an observable GW signal. It neither enforces the realized cosmological path nor computes a spectrum. The current runner invokes **MinimaTracer and CalcTemps, not CalcGW**. GW strength needs the transition dynamics, released energy and time scale, with assumptions about expansion and bubble-wall/fluid behavior; the field-space jump alone does not fix these quantities. [BSMPT v3](https://arxiv.org/abs/2404.19037).

### BSMPT setup and optional prefilter

The installed source is BSMPT **3.2.0 with the custom TRSM model**. The wrapper defaults to `Tlow=0` via BSMPT, `Thigh=300 GeV`, `multistepmode=default`, and multithreading enabled. The source has `C_UseParwani=false` (Arnold–Espinosa). The TRSM renormalization scale is the SM VEV. BSMPT uses its loop-corrected Landau-gauge potential; no gauge/scale/resummation uncertainty band is generated here. [BSMPT program description](https://arxiv.org/html/2404.19037v2#S3.SS2).

CalcTemps' NLO stability check is on by default. The high-temperature restoration check is on as a status check; the wrapper does not request an EWSR exclusion mode. A trace only covering up to 300 GeV does not establish a complete higher-temperature history. A plotted phase ending is not, by itself, a demonstrated first-order jump; coexisting minima and their potential crossing must be identified.

The optional `--ewpt-require-eq418` tests positive diagonal entries, all pairwise principal minors, and a determinant of a matrix built with `c_i=lambda_i/2`, `c_ij=lambda_ij/2`. This is an additional positive-definiteness-style prefilter, not the general copositivity condition and not a necessary FOPT criterion. The repository does not identify the original Eq. 4.18 reference sufficiently to authenticate its intended normalization. It is off by default; I would not use it as a general candidate-selection requirement without resolving that reference and purpose.

### Your plotted point 223615

| Diagnostic | Result from the archived payload |
|---|---|
| Critical EW-entry temperature | 154.147 GeV |
| Critical EW-entry jump/T | **0.13598** → baryogenesis flag **False** |
| Largest total-field critical jump/T | **6.59080**, at **51.5525 GeV** → GW flag **True** |
| Nucleation/percolation for that low-temperature transition | `not_met` |
| Nominal micrOMEGAs freeze-out temperature | about **0.7535 GeV** |
| Equilibrium Z2 compatibility | **True**, with the interpretation limits above |
| Maximum fractional EW/S VEV change near Tf | `2.26e-6` / `1.45e-6` |

Thus it is exactly the kind of **large singlet-driven transition candidate** your collaborator wants to retain. The coexistence of branches in the original plot is consistent with distinct local minima. Its equilibrium envelope favors the required X-unbroken vacuum at low temperature, but the saved transition history ends before the final transition, so the payload does not establish that the universe reached that final vacuum. It remains selected as a GW candidate under your agreed policy.

## 6. Validation, provenance and recommended order of work

### What was checked

- Read the active generator, theory and RG formulas, HiggsTools adapter, STU/W-mass code, DM model/driver/wrapper, CMB and DD adapters, BSMPT model/wrapper/status flow, candidate functions, thermal history functions, output schema and plot classifications.
- Inspected installed model defaults, HiggsTools version/datasets and BSMPT/micrOMEGAs source state on manto.
- Ran the counterexamples in [reproduce_constraints.py](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/audits/2026-09-23/reproduce_constraints.py). Results are in [numerical-probes.json](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/audits/2026-09-23/numerical-probes.json). The script's stdout may also contain native HiggsTools warnings, preserved in the accompanying `.log`.
- Evaluated the current candidate/history functions on the existing point-223615 payload; no fresh BSMPT run was needed for this demonstration.
- Ran **242 regression tests: 240 passed, one failed, one skipped**. The failure is `test_trsm_cmb.TestCMBPhysics.test_disabled_does_not_add_columns_or_change_results`: its assertion expects CMB columns at the very end, while the current writer appends the new EWPT candidate columns after them. Header and row construction use the same order; this is an out-of-date ordering assertion, not demonstrated numerical column corruption. The skipped test requires an archived generated fixture absent on manto. [Test log](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/audits/2026-09-23/regression-tests.log).

Passing regression tests validate the tested software behavior; they did not catch the scientific/convention gaps above. No corrected relic-density campaign, new tunnelling calculation, or complete fit/likelihood validation was performed in this audit.

### Installed provenance

| Component on manto | Audited state |
|---|---|
| HiggsTools | 1.1.3 |
| HiggsBounds dataset | `3539f9adf3724f8752701208bf059d1658824cd2`, commit dated 2025-09-11 |
| HiggsSignals dataset | `3c6c2538d9ba0549ea0559cf9b38c18eb80c3af1`, commit dated 2023-05-22 |
| BSMPT | `3.2.0-dirty`; model-registration changes; TRSM source matches the repository copy |
| micrOMEGAs | Default 6.1.15; optional 7.1.4 |

These identify the installed datasets, not a claim that they include all results available in September 2026. Scan metadata saves options, DM backend and DD-table provenance, but does not fully fingerprint repository code, HiggsTools datasets or the modified BSMPT build. Updating those components can therefore change resumed-scan physics without an immutable-options mismatch. Save their versions/hashes and the criterion definitions for future campaigns. [Runtime inventory](/Users/apapaefs/Projects/TwoSingletDM/twosingletDM/audits/2026-09-23/runtime-inventory.json).

### Recommended repair order

1. Align SM masses/VEVs and validate the DM photon/gluon effective vertices; quantify the change at resonance, cancellation and threshold benchmarks.
2. Add explicit zero-temperature vacuum assessment; define and implement the intended RG validity criterion; fix the mass initialization and use a symmetric unitarity eigensolver.
3. Separate failed/unassessed BSMPT results from negative candidates and stabilize the reconstructed global branch near crossings.
4. Retain your agreed critical-only baryogenesis and any-field GW cuts, with nucleation/percolation diagnostic rather than veto status. Rename the legacy “strong FOPT” plots.
5. Replace or clearly qualify the inherited DD fit, align the HiggsSignals reference, improve numerical precision/error reporting, and freeze provenance before reevaluating affected scans.

There is no scientific reason in this audit to impose a new nucleation/percolation hard cut on the qualitative candidate sample you requested. The necessary improvement is to make the physical meaning and calculation status of each saved flag explicit.
