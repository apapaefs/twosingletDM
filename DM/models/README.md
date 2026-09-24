micrOMEGAs uses the generated CalcHEP tables in `h4GOn/` (the production choice)
or `h4GOff/`. Their authored source is `lanhep_mdl/TRSM_mixed.mdl`.

Use [the setup script](../README.md) to install the model, `main.c` and
`trsm_loop.c` together in a fresh micrOMEGAs release. Copying the tables alone
is insufficient because their loop coefficients call the helper library.

[The v2 change audit](changes-v2.md) explains the physics changes and the
restoration of the historical two-Higgs operators. Both variants include
`h_i h_j gamma gamma` and `h_i h_j gg`; `h4GOn` also includes the restored
`x1/X1`, `G2t/G2T` auxiliary chain for `h_i h_j gggg`, alongside the retained
single-Higgs `G1t/G1T` chain. The installer selects `h4GOn`; no extra option is
needed to enable the restored vertices.

The restored contacts retain their historical `LAAh`, `LGGh`, `RQCDh`
coefficients. The corrected single-Higgs coefficients and annihilation hook
remain separate. This restores the previous effective approximation; it does
not provide a new matched two-Higgs loop calculation. Rebuild in a fresh
installation to use these tables; existing compiled installations are unchanged.

To regenerate both variants, from the repository root:

```bash
sh DM/tools/regenerate_models.sh /absolute/path/to/LanHEP/lhep
```

Review the generated diff and validate before rebuilding an installation.

The operator regression is `python -m unittest test_trsm_dm_model_normalization`.
For a native check, use a fresh, inactive installation built with the restored
tables, and run from the repository root:

```bash
TRSM_MODEL_DIR=/absolute/path/to/micromegas_7.1.4/TRSM
cp DM/tools/check_restored_operators.c "$TRSM_MODEL_DIR/"
make -C "$TRSM_MODEL_DIR" main=check_restored_operators.c
TRSM_RUNTIME_DIR=$(mktemp -d) \
  "$TRSM_MODEL_DIR/check_restored_operators" "$PWD/DM/data.par"
```

The same check supports 6.1.15. It checks all three Higgs-pair photon/gluon
contacts against their analytic effective-operator normalization, and exercises
the `x1/G2` chain with physical exchange diagrams excluded. It checks mixing
factors and independence from `Maux`, at a fixed 2 TeV phase-space point,
`SinT=0.3` and `Mh2=200 GeV`. This isolates the restored approximation; it is
not a prediction of a physical four-gluon rate or a full loop matching test.
