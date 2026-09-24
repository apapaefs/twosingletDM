# Constraint pipeline v2

## Scope and selection

`trsm_constraints_v2`, schema 2, covers the active **vx = 0** model. Every evaluated input has a row, including exclusions and unavailable results. Boolean verdicts are `True`, `False` or unavailable (`None` in Python, `null` in JSON, `nan` in TSV). Unavailable calculations must never become negative candidate verdicts merely because their transition list is empty.

The zero-temperature subsets are independent:

- `experimental_subset`: HiggsBounds, HiggsSignals and STU, plus the W-mass comparison where covered.
- `dm_subset`: the configured relic/DD/line/CMB prescription, without a lower relic-abundance cut.
- `vacuum_tree_global`: verified desired tree vacuum; degenerate global minima pass with an explicit status; unresolved flat families remain unavailable.
- `rg_integration_success`, `rg_bfb`, `rg_unitarity`: numerical completion and the two physical running tests are separate.
- `theory_strict_subset`: verified tree global vacuum and all RG checks passing through 1 TeV.

`ewpt_eligible = (thc is True) and (experimental_subset is not False)`. Neither DM, vacuum nor RG veto EWPT exploration. Incomplete experimental assessments remain eligible and are shown separately. `evo` survives only as a labelled legacy numerical-integration result. The main TSV is a complete ledger; campaign `viable_count` describes the tree + assessed experimental + DM subset, not the number of retained rows.

### Qualitative transition candidates

| Flag | Positive criterion |
|---|---|
| `ewpt_baryo_candidate` | EW-symmetric → EW-broken FOPT with ΔvEW/Tc **> 1** |
| `ewpt_gw_candidate` | Any FOPT with √(Δw1² + Δwx² + Δws²)/T **> 1**, at an available Tc, Tn or Tp |

EW entry uses `|w1,false| < 5 GeV <= |w1,true|`. Fields are mapped to symmetry-equivalent positive representatives before calculating jumps and matching minima; opposite-sign copies of one phase do not produce a spurious FOPT. Equality at one fails the strict cut. A critical baryogenesis candidate cannot be rescued by a nucleation-only enhancement. A GW candidate can be selected by any of the three temperatures.

Positive evidence is retained even when another calculation is incomplete. A negative baryogenesis verdict requires sufficient critical-temperature assessment. A negative GW verdict additionally needs assessed nucleation/percolation outcomes, since they could provide a positive result. `not_met` is an assessed negative temperature outcome; failed/missing calculations remain unavailable. Percolation and completion are diagnostics, never vetoes of a positive candidate.

These are qualitative screens of CalcTemps FOPT pairs. They are not sphaleron-rate calculations, predictions of baryon production or GW spectra, or proof that the Universe traverses a particular transition. The perturbative finite-temperature potential and Δv/T retain their gauge and approximation dependence.

## Shared inputs and DM calculation

The single source is `config/sm-inputs-v2.json`: M1 = 125.09 GeV, GF = 1.1663787e-5 GeV^-2, MW = 80.385 GeV, MZ = 91.1876 GeV. The pipeline derives v = 246.21965079413738 GeV, sin(thetaW) from the on-shell mass ratio, and e = 2 MW sin(thetaW)/v. Cards, default LanHEP parameters, BSMPT inputs, conversions and resonance diagnostics use this convention. Cards and numerical output retain double precision.

The default is **micrOMEGAs 7.1.4**; `--micromegas-version 6.1.15` selects the corrected comparison backend. Both use the same corrected model and cross-section hook. `TRSM_LEGACY_VIRTUAL_OFF=1` explicitly reproduces the previous virtual-W/Z setting; the default enables both. Effective settings, pole masses, running-mass inputs, widths, solver status, Xf and Tf are saved. Failed relic calculations have no usable Xf/Tf or rescaled verdict; raw solver evidence and the reason are retained.

`DM/trsm_loop.c` uses complex one-loop fermion/W amplitudes at the actual sqrt(s), with coherent h1/h2 interference. Each mediator has one doublet projection. On-shell absolute coefficients give correctly normalized γγ and gg decay widths; the associated single-Higgs three/four-gluon operators use the same coefficient. The calculation uses effective running light/charm/bottom/strange quark inputs, pole top and W masses, and alpha_s at the loop energy. This single-Higgs prescription is leading loop, without the former heavy-top K factor. The independent validation includes both raw CalcHEP γγ/gg normalizations at a pure-mediator pole.

The historical two-Higgs γγ/gg contacts and the `x1/X1`, `G2t/G2T` four-gluon auxiliary chain have been restored separately, with their original `LAAh`, `LGGh`, `RQCDh` coefficient prescription. They retain the previous effective approximation; the single-Higgs validation does not establish a matched two-Higgs loop calculation. See the [model audit and restoration notes](../DM/models/changes-v2.md). Use a fresh build for this model revision; existing compiled installations and historical validation receipts are not updated automatically.

Every worker gets an isolated writable CalcHEP runtime/cache. Existing campaign installations remain intact. The driver refuses execution without the isolated runtime directory; the Python wrapper creates it automatically.

Approximate resonance sampling means **2 M3 ≈ M1 or 2 M3 ≈ M2**. The separate `--mass-ratio-m3-2m2 --delta-res ...` mode is a mass-ratio study. Metadata records the actual supported sampling envelopes. Exact poles are not rounded to four significant figures.

## Experimental prescriptions

- The default observed SI table is the [published LZ WS2024 release](https://doi.org/10.17182/hepdata.155182.v2/t1), with original and normalized data in `DM/data/lz2024`. It is a power-constrained observed 90% CL, per-nucleon isoscalar limit. Use log-log interpolation over 9–10,000 GeV; no extrapolation. DD outside coverage is unassessed. The explicit 2026 approximate-table option is retained.
- The relic upper bound is **0.121**, with abundance reference **0.12**. Rescaling is ξ = min(1, Ωh²/0.12), linear for DD and squared for gamma lines and CMB. The inherited Fermi R16 line prescription remains; its coverage/availability is recorded. CMB is enabled by default with version 7 and uses the backend's low-velocity s-wave prescription.
- HiggsSignals uses the **same widths/BR pipeline for the SM reference and scan points**. The rule remains Δχ² < 4 relative to that reference. It is a fixed screening prescription, not a profiled multi-parameter 95% confidence region.
- The correlated STU fit remains the Jens Erler private-communication input dated 16 May 2025, with χ² <= 7.82. Its oblique-parameter normalization retains alpha(0) = 7.2973525693e-3. The input fit is not authenticated as a public likelihood.
- The Tania Robens Snowmass W-mass curve is used only for 133–999 GeV, with its existing cubic interpolation. Coverage is explicit; STU and this curve are not presented as independent likelihood factors. Removable EWPO singularities, including M2 = MZ, are implemented; nonfinite results are unavailable.
- Physical Γ/M and cτ in mm are recorded independently of HiggsTools' provider width floor. No new broad-width or prompt-decay acceptance cut is implied.

## Vacuum and RG diagnostics

The canonical potential is evaluated on all stationary supports: origin, three axes, three planes and the full branch. Physical nonnegative squared fields are compared by potential depth and Hessian eigenvalues. Singular systems are assessed for feasible flat families and unbounded flat directions; unresolved cases are explicit. This is a tree assessment. BSMPT's NLO stability status is stored separately.

Unitarity uses the real symmetric scattering matrix and `eigvalsh`, avoiding complex artifacts from repeated polynomial roots. The one-loop running starts at 91 GeV with the inherited running SM couplings and tree scalar boundary values. The missing M3² term in the X mass parameter is repaired. BFB and unitarity are tested along the solution through **1000 GeV**, with the first detected violation bracket refined to 1e-5 GeV. Integration uses rtol = 1e-8, atol = 1e-10, solver points plus 513 logarithmic samples and a maximum step of (1000−91)/256 GeV. This is a numerical trajectory screen with **no threshold matching or precision pole-to-running matching**.

## Thermal interpretation

Executable return codes, NLO stability, tracing, coexistence and temperature statuses are separate, following the [BSMPT status documentation](https://arxiv.org/html/2404.19037v2#S3.SS12). The raw CalcTemps row preserves its calculated transition-history output. A failed independent MinimaTracer/refinement step does not erase valid CalcTemps evidence.

`PhaseProbe` refines competing local minima and evaluates their potentials at common temperatures. Symmetry copies are merged. A default crossing bracket is at most **1e-3 GeV**, when numerical accuracy permits; unresolved ordering and brackets remain explicit. “Resolved equilibrium” means ordering **among traced basins**, not a proof that every minimum was discovered. It is kept distinct from a cosmological history, especially if a transition never nucleates.

The separate freeze-out diagnostics are:

1. `dm_relic_z2_freezeout_compatible`: necessary equilibrium Z2 compatibility around/after nominal freeze-out. Unknown low-temperature ordering cannot pass this flag. It does not prove cosmological transition completion or validate the full relic calculation.
2. Thermal VEV variation in **[Tf/2, 2Tf]**, with the existing 10% screen and EW/X/S field thresholds 5/1/1 GeV.
3. Fixed-mixing proxies on compatible wx = 0 branches:
   K133(T) = [λPhiX w1(T) cos(a12) − λSX ws(T) sin(a12)]/2,
   K233(T) = [λPhiX w1(T) sin(a12) + λSX ws(T) cos(a12)]/2.
   Absolute changes, fractional changes and cancellation indicators are separate. Fractional changes are unavailable for a vanishing/cancellation-dominated reference denominator. These are sensitivity proxies; they do not compute thermal eigenstates, thermal masses or an updated Boltzmann evolution.
4. |Mi − 2 M3|/Tf and |Mi − 2 M3|/Γi, using the backend mediator widths. Near-resonance flags are diagnostics, not new exclusions.

## Outputs, provenance and plots

Generation, reevaluation, EWPT reprocessing, campaign aggregation and plotting use the shared v2 schema. Legacy fields remain readable and explicitly labelled. Old data cannot acquire v2 candidate flags without reevaluation. Every evaluated row retains independent verdicts, values and reasons. `ewpt_transition_strengths` preserves same-transition temperature pairs for plotting.

Manifests fingerprint the source commit and active source files, generated model, executable hashes, SM inputs, running prescription, Higgs datasets, widths/coupling tables, DD data and numerical settings. Incompatible resumes are rejected before appending output. Campaign aggregation reads actual output manifests rather than guessing filenames; `candidate_points.tsv` and `candidate_counts.json` use the new flags. `best_points.tsv` remains explicitly labelled as the legacy EW-jump ranking.

The plot suite adds parameter maps with DM/vacuum/RG overlays, same-transition Tc-versus-Tn/Tp jumps with missing counts, local-minimum trajectories and common-temperature potential differences, Tf/restoration brackets, and VEV/coupling sensitivity versus width-normalized resonance gaps. Plots display unavailable assessments separately.

## Build and run

Use a fresh `runtime-v2` sibling of the repository (or set `TRSM_RUNTIME_ROOT`). Preserve older installations. The runtime receipt is `runtime-manifest.json` in that directory.

1. Create a Python 3.13 virtual environment; install `config/runtime-requirements-v2.txt`. Build HiggsTools from the recorded source and use the recorded HB/HS dataset revisions. On these macOS hosts the local build used `CMAKE_ARGS='-DCMAKE_CXX_FLAGS=-D_LIBCPP_TEMPLATE_VIS='`.
2. Extract the original micrOMEGAs archives separately. The 7.1.4 archive SHA256 is `c8cf207b17541a5b7d7e7ff157f1c1eb36e1dd3f7547e7745559b195c2a34264`. Run `sh DM/setup_micromegas.sh <fresh-install>` for each version. The script refuses an existing TRSM installation. Regenerate model tables with `sh DM/tools/regenerate_models.sh <LanHEP/lhep>` when the LanHEP source changes.
3. Use BSMPT **3.2.1**, with the TRSM implementation at upstream source revision `32ccb11e856631f770e7e7471a85f4928729510b`. Run `python tools/stage_bsmpt_v2.py <upstream-BSMPT> <runtime-v2/BSMPT-3.2.1>`. Configure with the site's Conan toolchain and `-DCMAKE_BUILD_TYPE=Release -DBSMPTBuildExecutables=ON -DBUILD_TESTING=OFF`; build `CalcTemps MinimaTracer PhaseProbe Test`. The staging script checks the version; the runtime receipt fingerprints the full compiled source tree, including the base BSMPT kernels.
4. Run `tools/runtime_manifest.py` with the BSMPT and HiggsTools source locations. It verifies installed model/driver/helper files before writing the receipt.
5. Use `tools/run_next_scan.py` with `config/next-scan-v2.json`. It prints the exact command; `--run` launches in a fresh campaign directory. `--pilot --run` limits generation to four draws across two isolated workers. Resume individual seed outputs with `generate_trsm_points.py --resume-from <actual-main-path> --nrandom <total-target>`.

The BSMPT TRSM source revision is based on public upstream `431df3b6a39cd5458e70f0035f6be87c40d3394d` (temperature-derivative and sound-speed fixes). The model overlay needed on top of that source is tracked under `BSMPT/` in this repository. Manto also retains the canonical source checkout at `runtime-v2/BSMPT-upstream-32ccb11e`.

The supplied next-scan configuration uses 10,000 draws across four seeds, a 2 GeV annihilation-pole window and Tmax = 1000 GeV. Existing mass, VEV, mixing and portal ranges are written explicitly. This is a reproducible starting configuration, not an optimization of scan coverage. The full scan is not launched by validation.
