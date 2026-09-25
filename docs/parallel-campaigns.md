# Parallel scan campaigns

The campaign launcher runs independent Python scans on one machine. Each scan
has a distinct seed, private writable files and its own checkpoint. The number
of scans, stopping target and maximum number running at once are separate
settings. The full v2 workflow runs in each worker, including EWPT for eligible
points. All evaluated points and the independent constraint flags are retained.

## Installation and launch

Clone the repository and follow [the v2 runtime setup](constraints-v2.md) to
compile/install micrOMEGAs, BSMPT (including PhaseProbe), HiggsTools and the
HiggsBounds/HiggsSignals datasets. Set `TRSM_RUNTIME_ROOT` if the compiled runtime
is outside the default sibling `runtime-v2` directory. The datasets retain the
documented sibling `hbdataset`/`hsdataset` layout. Activate the Python environment
containing the scientific dependencies before using the commands below.

For example, from the repository:

```bash
python tools/run_next_scan.py \
  --campaign-dir output/my-campaign \
  --seed-start 10000 --nseeds 200 --nrandom 5000 \
  --nrandom-count-evo-thc --jobs 192 \
  --run
```

These numbers are examples. Omit `--run` to print and validate the command without
launching anything. `--jobs auto` uses the CPUs available to the process. Explicit
job limits are also capped by available CPUs and the number of remaining scans.
On Linux, CPU affinity is respected. Choose a lower limit when RAM or other work
on the host requires it. CPU detection does not estimate a safe RAM allocation.

`--nrandom-count-evo-thc` means **each seed stops after `--nrandom` points with
both `evo=True` and `thc=True`**. It does not require the newer global-vacuum or RG
subset flags, DM acceptance, or an EWPT candidate. Those flags remain available
for later filtering. More raw draws may be needed to reach the target, and all
those rows remain in the output. Without this option, the target counts raw
draws. A completed scan means its target was reached, not that every physics
calculation returned an assessed result.

The wrapper reads `config/next-scan-v2.json`. CLI values override configuration.
The existing `draws_per_seed` key remains supported; new configurations can use
`points_per_seed` and the Boolean `nrandom_count_evo_thc`. Other keys include
`seed_start`, `nseeds`, `jobs` (a positive integer or `"auto"`),
`checkpoint_every`, `heartbeat_seconds`, `shutdown_grace_seconds`,
`ewpt_thigh_GeV`, and `generator_arguments`. The current supplied profile and its
raw-draw defaults are unchanged. The wrapper always enables EWPT; the lower-level
`run_trsm_seed_campaign.py` additionally supports campaigns without EWPT.

A fresh campaign validates the selected Python interpreter, HiggsTools datasets,
micrOMEGAs v2 capabilities and requested BSMPT executables before launching its
workers. `preflight.log` contains diagnostics. Explicit runtime paths can still
be provided through the existing generator arguments; the wrapper also exposes
`--python-executable`, `--ewpt-executable` and `--ewpt-minima-executable`.

Each worker uses `--no-ewpt-multithreading`, single-thread numerical-library
settings, and single-thread HiGHS solves. BSMPT's C++ thread pool must be disabled
explicitly: setting `OMP_NUM_THREADS` alone does not control it. Standalone
`generate_trsm_points.py` retains its default BSMPT multithreading and accepts
`--ewpt-multithreading` / `--no-ewpt-multithreading` to select it.

## MadGraph cross sections

Both launchers accept `--run-mg5`. By default this evaluates `gg_heta0` and
`pp_eta0Z` at 13.6 TeV for points passing full viability, including flavour.
Repeat `--mg5-process NAME` to choose from `hh`, `hhh`, `gg_heta0`, and `pp_eta0Z`.
`--mg5-without-dm` drops only the DM requirement; all non-DM constraints,
including flavour, still apply. These options do not change EWPT eligibility.

For ten scans of 100 evo/thc-passing points using the sampling of the historical
Manto seed-17999 scan:

```bash
source ../runtime-v2-new/activate.sh
python run_trsm_seed_campaign.py \
  --campaign-dir output/odysseus-smoke-10x100-mg5 \
  --seed-start 17999 --nseeds 10 --nrandom 100 \
  --nrandom-count-evo-thc --jobs 10 \
  --run-ewpt --ewpt-thigh 1000 --run-mg5 \
  --generator-extra-arg=--independent-m3 \
  --generator-extra-arg=--scan-k133-k233-log \
  --generator-extra-arg=--m3-max=65 \
  --generator-extra-arg=--no-print-info
```

The current v2 physics and flavour prescriptions apply. All raw draws remain in
the ledger; the target counts evo/thc passes. A scan can finish without finding
an MG5-eligible point. Cross sections appear in `mg5_xsec_<process>_pb` columns
in each scan and `combined_points.tsv`.

Activate the runtime before launching so `TRSM_MG5_LOCATION` selects the compiled
processes. Preflight checks the selected launchers and compiled subprocesses and
records stable process/model files in the campaign receipt. Missing installations
stop the campaign before workers start. Actual run cards and helicity-optimized
binaries are mutable and are not included in this receipt's source hashes.

MG5 updates cards and event files inside each generated process directory.
An exclusive file lock covers each invocation: workers using the same process
take turns, while different processes and other scan calculations can proceed
concurrently. Campaign workers set `TRSM_MG5_CORES=1`, so each invocation uses one
core. `--jobs` controls concurrent scans, not independent copies of MG5.

The `tools/run_next_scan.py` wrapper exposes the same options. Its configuration
also accepts `run_mg5` and `mg5_without_dm` booleans and an `mg5_processes` list;
CLI values override them. MG5 is disabled by default. Supply MG5 options directly
to the launcher, not through `--generator-extra-arg`.

MG5 enablement, its DM selection, and its process list are saved with the
campaign and cannot change on resume. Use a fresh campaign directory to enable
MG5 for a previous scan that did not run it. Existing native installations can
be reused; this feature does not require a rebuild.

## Files, progress and interruption

```text
campaign_state.json                 effective configuration, receipt, seed ledger
preflight.log                       installation validation
logs/seed_<seed>.log                 appended attempt logs
manifests/seed_<seed>.json           actual generator output paths
seeds/seed_<seed>/output/            TSV, metadata and checkpoint
seeds/seed_<seed>/cache/             private plotting cache
ewpt/seed_<seed>/point_<index>/      retained EWPT results
campaign_summary.json/.tsv          seed states and raw/evo-thc/viable counts
combined_points.tsv                 all checkpoint-committed evaluated rows
candidate_points.tsv                baryogenesis or GW candidate rows
candidate_counts.json               separate true/false/unassessed counts
best_points.tsv                     historical EW-jump ranking
```

`--run-cwd` on the lower-level launcher optionally selects a different **parent**
for private `seed_<seed>` working directories. It no longer shares one working
directory between seeds. Phase-plot basenames must be relative to each point's
EWPT directory. The DM/EWPT installations and datasets are shared read-only;
micrOMEGAs/CalcHEP temporary workspaces remain private to each worker. MG5 uses
the locked, writable process directories described above.

Heartbeats distinguish queued, running, complete, failed and interrupted scans,
and report raw draws and evo/thc passes from checkpoints. Use
`--heartbeat-seconds` to change the interval; zero disables heartbeats.

Ctrl-C or SIGTERM stops new launches and signals the campaign's own worker
process groups, including native children. Workers have 30 seconds by default to
exit and checkpoint before remaining owned processes are killed. Change this
with `--shutdown-grace-seconds`. `--checkpoint-every N` controls checkpoint
cadence, defaulting to one completed draw. Recovery may redo uncommitted draws.
Never remove locks or checkpoint files to work around a live-worker error.

## Resume and aggregation

```bash
python tools/run_next_scan.py --campaign-dir output/my-campaign \
  --resume --jobs 96 --run

python tools/run_next_scan.py --campaign-dir output/my-campaign \
  --aggregate-only --run
```

Resume loads the saved effective configuration, skips complete seeds, continues
unfinished seeds at their saved RNG states, and starts queued seeds. Logs are
appended. Worker failures do not cancel other seeds. An incomplete campaign
returns a nonzero status and lists its incomplete seeds in the JSON summary.

The concurrency, heartbeat interval, checkpoint cadence and shutdown grace can
change on resume; omitted values retain the saved execution settings. Seed
identities, stopping target/counting mode and physics configuration cannot
change. The source/runtime receipt is checked before resuming. Use the original
source revision and installation for existing runs; a changed physics source or
runtime requires a fresh campaign. Cross-host relocation of existing checkpoints
is outside this workflow: clone and compile on Odysseus, then start a fresh
campaign there.

A supervisor lock and per-scan locks reject duplicate writers. If a supervisor
was killed while workers survived, wait for those workers to finish before
resuming. Recovery uses actual manifests or the unique checkpoint in a seed's
private directory, so midnight and missing manifests do not cause fresh runs to
replace existing output. Output lacking a usable checkpoint is reported for
inspection rather than silently overwritten.

Aggregation needs no scientific runtime. It streams data in seed/point order,
checks each committed prefix, excludes unfinished output tails and uncommitted
EWPT directories, and identifies incomplete seeds. It does not modify the scan
TSVs. Repeating aggregation rebuilds the products without duplicating rows.
`--aggregate-only` also takes the supervisor lock; run it after stopping the
campaign. The Python API's `CampaignResult.combined_points` is a re-iterable TSV
view with a row count, rather than an in-memory list.

## Validation and sizing on Odysseus

Run the fast regression tests:

```bash
python -m unittest test_parallel_campaign test_run_trsm_seed_campaign \
  test_trsm_scan_resume test_generate_trsm_points_runner test_trsm_ewpt_runner \
  test_constraint_v2 test_trsm_DM_unit test_trsm_theory_constraints_unit
```

Run the small native pilot after compilation:

```bash
python tools/pilot_parallel_scan.py --output-dir output/parallel-pilot
```

It runs the same two seeds serially and concurrently, targeting one evo/thc pass
each in a known EWPT-eligible smoke region. It compares sampled coordinates,
constraint flags and numerical results using the existing DM/Higgs/EWPT
validation tolerances; eligible points must complete both MinimaTracer and
CalcTemps. It samples native thread counts, worker-process-group RSS, elapsed
time, throughput and directory size into `report.json`. `ps` must be permitted
to inspect the pilot's own processes. Summed RSS includes native children and
can count shared pages more than once; it is a sampled measurement, not an
exact memory reservation.

`--nseeds`, `--nrandom`, `--jobs`, `--seed-start`, `--sample-seconds` and
`--timeout-seconds` are configurable. Supply `--config` to measure the intended
production ranges after the smoke pilot; that sample must include EWPT-eligible
points to pass. Increase concurrency gradually while checking memory and disk
usage. The small smoke region establishes correct execution, not representative
production throughput or an optimal worker count.
