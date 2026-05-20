# BSMPT TRSM integration

This directory archives the BSMPT files generated for the TRSM implementation with
the input basis

```text
m1 m2 m3 vs a12 lx lphix lsx
```

where `m1`, `m2`, and `m3` are the scalar masses, `vs` is the singlet vev,
`a12` is the doublet-singlet mixing angle, and `X` has zero vev.

The directory mirrors the relevant BSMPT paths:

- `include/BSMPT/models/ClassPotentialTRSM.h`
- `src/models/ClassPotentialTRSM.cpp`
- `tools/ModelGeneration/Mathematica/TRSM.wl`
- `example/TRSM_Input.tsv`
- registration files under `include/BSMPT/utility/` and `src/models/`

Validation was run in the sibling BSMPT checkout with:

```bash
cmake --build build/macos-armv8-release -j 8
build/macos-armv8-release/bin/Test --model=TRSM --input=/Users/apapaefs/Projects/TwoSingletDM/BSMPT/example/TRSM_Input.tsv --line=2
```

The model test reported `21 tests out of 21 passed`.
