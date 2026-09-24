# Dark matter implementation

The active `vx=0` model uses the corrected v2 driver (`main.c`), coherent loop
correction (`trsm_loop.c`, `trsm_loop.h`) and generated model in `models/h4GOn/`.
Use a rebuilt installation containing all of these components; copying only
`main.c` into an older installation is insufficient.

- [Standalone examples and reproducible CMB benchmarks](steer_example/README.md)
- [Planck CMB prescription](planck-cmb.md)
- [Direct-detection table and conventions](direct-detection.md)
- [Constraint profile and runtime setup](../docs/constraints-v2.md)
- [Validated release and scan commands](../docs/next-scan-readiness-v2.md)

The default backend is micrOMEGAs 7.1.4. Corrected 6.1.15 is retained for explicit
comparisons. The standalone example reads the same limits and comparison code
as the scan and preserves raw output for inspection.
