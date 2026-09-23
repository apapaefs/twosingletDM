# Constraint v2 implementation record

The approved active-vx=0 pipeline is implemented and validated. The canonical
specification and readiness evidence are in `docs/constraints-v2.md` and
`docs/next-scan-readiness-v2.md`.

## Completed

- Shared full-precision SM inputs, regenerated LanHEP models, corrected coherent
  AA/GG loop hook, micrOMEGAs 7.1.4 default and corrected 6.1.15 comparison.
- Independent tree-vacuum and running subsets; dense RG trajectory assessment
  through 1 TeV; stable scattering-matrix eigenvalues.
- Verified observed LZ WS2024 numerical release, unified abundance rescaling,
  common HiggsSignals reference, EWPO limits and physical-width diagnostics.
- Nullable EWPT candidates, common-temperature minimum refinement, separate
  equilibrium/cosmological interpretation and thermal relic-input diagnostics.
- Complete point ledger, provenance/resume protection, actual campaign manifests,
  candidate summaries and all five requested plot families.
- 343 passing tests on each host, independent amplitudes, 28 runtime cases,
  audit counterexamples, deterministic 100-point comparison and bounded pilots.

## Runtime consistency

The final provenance check detected different BSMPT base kernels on the two
hosts. Both release installations now use the canonical BSMPT 3.2.1 source;
point 223615 and the integrated pilots were repeated. Active pipeline sources,
BSMPT/HiggsTools kernels, model sources, SM inputs and HB/HS datasets match.
All 100 zero-temperature numerical comparisons match at the saved precision.

Final runtime receipts are generated against the committed source after
production deployment. Existing campaign installations and outputs remain
preserved. The 10,000-point next scan has not been launched.
