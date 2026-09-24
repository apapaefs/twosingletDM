# Audit of the micrOMEGAs model changes

The changes questioned here are in commit **`34beddb`** (23 September 2026),
before the parallel-campaign work. The standalone example rewrite is the
separate commit **`33e39f1`**. The initial interface compatibility repair left
the physics files unchanged. Following the explicit request to reinstate the
removed operators, the historical two-Higgs contacts and auxiliary chain are
now restored as described below.

## `x1` and the four-gluon interaction

Auxiliary fields are required to express these contact interactions in
CalcHEP's vertex format. The collaborator's concern about keeping that
representation is valid. However, the old model uses two distinct auxiliary
chains. Inspection of the actual old vertices gives:

| Chain | External interaction represented |
|---|---|
| `h_i G G G1t.t` joined to `G G G1T.t` | One Higgs plus four gluons |
| `h_i h_j x1`, `X1 G G G2t.t`, `G G G2T.t` | Two Higgs fields plus four gluons |

The first chain remained in `h4GOn/prtcls1.mdl` and `h4GOn/lgrng1.mdl`, for
both `h1` and `h2`. The second was removed along with the two-Higgs effective
loop operators and is **now restored**. Thus removing this particular `x1`
did not remove the single-Higgs four-gluon operator. `h4GOff` omits the
auxiliary three/four-gluon representation; the installer selects `h4GOn`.

This distinction is also visible in the old LanHEP source. In unitary gauge,
with the neutral doublet fluctuation `h = CosT*h1 + SinT*h2`,

```text
(shd*shD - vevh^2/2)/vevh = h + h^2/(2*vevh)
(shd*shD - vevh^2/2 - vevh*h)/vevh = h^2/(2*vevh)
```

It is the **second**, purely quadratic expression that multiplies `x1*Maux`.
The auxiliary fields have constant propagators and are not new physical
particles; see the [LanHEP auxiliary-field documentation](https://theory.sinp.msu.ru/~semenov/doc/node10.html).

## What changed and why

| Change in `34beddb` | Reason and scope |
|---|---|
| Replace `LAAh=-abs(lAAhiggs(Mh,"h1"))` and `LGGh=-abs(lGGhiggs(Mh,"h1"))` with separate `LAA1/2`, `LGG1/2` | The old helper already contained the h1 doublet projection, which the vertices multiplied again; it also used the h1 loop scale for h2. The new coefficients carry one projection and use each mediator's mass for on-shell widths. |
| Add `DM/trsm_loop.c` and its `improveCrossSection` hook | For `XX -> gamma gamma, gg`, evaluate the complex one-loop form factor at the actual annihilation energy and add the h1/h2 amplitudes coherently. An absolute on-shell coefficient cannot describe their general off-shell interference. |
| Remove `RQCDh` | The implementation now states a leading-loop prescription rather than applying the previous heavy-top inclusive correction factor to all kinematics. This is a perturbative-prescription change. |
| Replace the doublet bilinear effective operator by terms linear in `h1/h2` | This removes the approximate `h_i h_j gamma gamma` and `h_i h_j gg` contact operators, and the associated `x1/G2` chain. A single-Higgs loop form factor does not by itself define the two-Higgs box/contact amplitudes. This is a **reduction in model scope**, not a required consequence of fixing the extra mixing factor. |
| Retain `G1t/G1T` and use `LGG1/2` in its vertices | Preserve the effective single-Higgs three/four-gluon completion. This is not a full finite-mass one-loop calculation of arbitrary multi-gluon processes. |
| Change `EE`, `SW`, `Mh`, `Mtp` defaults | Align the DM calculation with the shared v2 inputs: GF/MW/MZ electroweak scheme, Mh=125.09 GeV, top pole mass=172.5 GeV. Actual cards carry full precision. |
| Introduce `B00000...B00015` and rewrite some fermion-Z expressions | LanHEP regeneration abbreviates long scalar expressions and rearranges algebra. The canonical scalar/DM potential is unchanged by this regeneration. |

The earlier canonical X-potential normalization is a different change,
commit **`3d54f24`**. Its factors must not be conflated with the later removal
of the loop-generated two-Higgs operators.

The scope reduction should have been highlighted clearly when it was made.
That reduced model was not backward compatible for calculations using those
two-Higgs operators. The restoration includes their vertices and coefficient
prescription, as well as the auxiliary particle declarations.

## Requested restoration

The LanHEP source and both regenerated CalcHEP variants now restore all the
vertex rows removed by `34beddb`: six two-Higgs photon/gluon contacts in each
variant, plus the five `x1/X1`, `G2t/G2T` chain vertices in `h4GOn`. Their
coefficients and Lorentz structures match the historical tables exactly.
All vertex rows that were already present before this restoration are
unchanged. No additional run option is required.

To avoid introducing another physics prescription while restoring compatibility,
the reinstated terms use the original definitions:

```text
LAAh  = -abs(lAAhiggs(Mh, "h1"))
LGGh  = -abs(lGGhiggs(Mh, "h1"))
aQCDh = alphaQCD(Mh)/pi
RQCDh = sqrt(1 + 149/12*aQCDh + 68.6482*aQCDh^2 - 212.447*aQCDh^3)
```

These multiply the purely quadratic expression
`(Phi†Phi - v^2/2 - v*h)/v = h^2/(2v)` and the associated auxiliary chain.
Their historical h1 form-factor scale, mixing dependence, signs and QCD
factor are preserved. They are the inherited effective approximation, not
a newly matched two-Higgs loop result. In particular, this does not introduce
a new two-Higgs/three-gluon vertex absent from the historical CalcHEP tables.

The linear `LAA1/2`, `LGG1/2` vertices and the complex `XX -> gamma gamma, gg`
annihilation hook retain the v2 prescription. Restoring the old coefficients
for the quadratic terms does not reintroduce them into the linear terms.
The X potential, SM inputs and scan settings are unchanged.

Regenerate with `sh DM/tools/regenerate_models.sh /path/to/LanHEP/lhep`.
Install with `sh DM/setup_micromegas.sh /path/to/fresh/micromegas_7.1.4`
(or `6.1.15`), then select that executable using `--micromegas-main`.
Existing compiled installations, caches and campaign outputs are left intact;
they do not acquire these vertices until a new model is built. Historical
validation receipts describe their original model revision.

## Evidence and limits of the validation

`test_trsm_dm_model_normalization.py` checks the scalar/DM normalization,
the retained single-Higgs chains, all restored contact and auxiliary vertices,
the `h4GOn`/`h4GOff` distinction and the separate coefficient prescriptions.
`tools/validate_dm_v2.py` compares
the compiled loop amplitudes with independent high-precision formulae for
both supported backends: pure/mixed mediators, poles, W/top thresholds,
interference cancellation, on-shell widths, and raw CalcHEP gamma-gamma/gg
normalization. Those 22 cases pass when rerun for this audit.

The restoration was built in fresh micrOMEGAs 6.1.15 and 7.1.4 installations.
`DM/tools/check_restored_operators.c` passes nine native process checks per
backend: all three Higgs pairs into photons, gluons, and the isolated four-gluon
auxiliary chain. The two-vector contacts agree with their analytic effective
vertex normalization. The chain has the expected Higgs-pair mixing dependence
and is independent of the arbitrary `Maux` parameter (tested at 1 and 7).
Build/run instructions are in [the model README](README.md).

The 35 focused DM/CMB/model tests and the four native 7.1.4 DM benchmarks
(86 checks at the existing tolerances) also pass. These checks establish
restoration of the historical approximation and preservation of the tested
single-mediator-loop DM prescription. They do not establish full two-Higgs
loop amplitudes or independently validate all multi-gluon differential
amplitudes. The isolated auxiliary contribution is not a physical cross section.

Inspect the original changes from the repository root:

```bash
git show 34beddb -- DM/models DM/trsm_loop.c DM/main.c
git show 34beddb^:DM/models/lanhep_mdl/TRSM_mixed.mdl
python -m unittest test_trsm_dm_model_normalization
python tools/validate_dm_v2.py validation/my-dm-amplitude-check
```

The native check needs both compiled backends and `mpmath`; it is optional for
ordinary scans. Keep historical and current installations separate when
comparing results, and retain their cards, executable identities and raw logs.
