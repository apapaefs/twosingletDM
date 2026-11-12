# MG5 survey repair validation — 25 September 2026

The reported zero was reproduced at smoke-campaign seed 18008, draw 93.
The old wrapper requested one iteration while the unsplit MG5 survey required
three; channel estimates were nonzero but final survey files contained zeros.
Raising the limit to three exposed a second failure in MG5's adaptive helicity
sampler. Explicit helicity summation (`nhel=0`) completed successfully.

## Frozen calculation

MG5_aMC 3.5.15 on Odysseus, GNU Fortran 15.2.0, existing
`loop_sm_twoscalar_generic` model, proton beams at 13.6 TeV and bundled `nn23lo1`
PDFs. The original failure was in loop-induced `g g > h eta0 [noborn=QCD]`.
M2 = 511.1 GeV, M3 = 61.66 GeV, Gamma2 = 2.29 GeV.
The full parameters, widths, scales and run settings are retained in
[point.json](native/point.json) and [original banner](native/original-banner.txt).
The repair changes integration settings, not the model or physics parameters.

The old and updated wrappers used separate process copies with identical inputs
and random seed 21. The independent reference used seed 22, 3000 requested survey
points and a maximum of six iterations, compared with 300 points and three
iterations in the corrected scan. All runs used one core and requested one
unweighted event. Each successful run had exactly one event, finite positive
normalization and event weight, consistent survey bounds, and exact loop
helicity sums. The additional hh, hhh and eta0Z checks are smoke tests at the
same parameter point, not independent precision validations of those processes.

## Native results

| Case | Process | Cross section [pb] | Integration error [pb] | Wall time [s] |
|---|---|---:|---:|---:|
| current | gg_heta0 | 6.68066e-05 | 8.067224e-07 | 6.62 |
| legacy-repaired | gg_heta0 | 6.68066e-05 | 8.067224e-07 | 6.57 |
| reference | gg_heta0 | 6.72717e-05 | 2.292885e-07 | 43.58 |
| current | pp_eta0Z | 7.7661e-05 | 9.1017e-07 | 2.52 |
| current | hh | 0.03117 | 0.0003516255 | 10.18 |
| current | hhh | 9.87392e-05 | 2.785908e-06 | 227.55 |

The corrected wrapper and repaired original wrapper gave identical cross
sections and integration errors. The independent reference differs by
0.555 combined integration standard deviations (required tolerance: 3).
The one-iteration false zero and the sampled-helicity failure are retained as
[baseline log](native/baseline-gg_heta0.log) and
[three-iteration sampled-helicity log](native/fixed-gg_heta0.log).
All successful logs and machine-readable results are in `native/`.

## Deployment and recovery

The three installed loop processes on Odysseus were repaired under their writer
locks while both saved campaigns were locked and confirmed stopped. Each repair
changed only `bin/internal/gen_ximprove.py`, `Cards/run_card.dat` and
`Cards/run_card_default.dat`. Originals and per-file SHA-256 receipts are kept
under each process's `.trsm-maintenance/survey-minimum-*/` directory. The guard
preserves legitimate one-iteration split-grid chunks.

Full preflight checks before and after the repair returned identical receipts,
matching both saved campaigns. All 3547 protected source, physics-input,
event, log and campaign-output files retained identical hashes. No checkpoint
or fingerprint was rewritten. Neither campaign was restarted; no pilot worker
was left running. Details: [repair report](native/repair-report.json).

The source changes are local and uncommitted/unpushed. The original Odysseus
checkout remains clean at `c29d1c28d4d39f84168499c8feaadc7a0cd75511`, preserving
resume compatibility. Resume from that original checkout and activated runtime:

```bash
python run_trsm_seed_campaign.py \
  --campaign-dir output/odysseus-100x1000-66661 --resume --jobs 100
```

The repair tool and generic instructions are documented in
[the runtime guide](../../../docs/runtime-bootstrap.md#repairing-an-existing-mg5-survey).
`pilot.py` and `deploy.py` are installation-specific audit scripts, not campaign
entry points. Physics-model changes and campaign restarts are outside this fix.

## Automated regression checks

108 tests passed across `test_madloop_build`, `test_bootstrap_runtime`,
`test_mg5_process_runner`, `test_mg5_survey_repair`, `test_parallel_campaign`,
`test_trsm_scan_resume` and `test_run_trsm_seed_campaign`.
The new checks cover survey limits, loop-only helicity configuration, preserving
zero-survey failures, backups and idempotence, active-writer rejection, unknown
source/card refusal and out-of-process source links. `git diff --check` passed.
