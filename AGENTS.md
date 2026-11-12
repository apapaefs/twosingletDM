# Project working agreements

These rules apply to ASTRA and other coding agents working in this repository.
Prioritize scientific intent, reproducibility, backwards compatibility and
readable documentation over architectural elegance. Preserve existing user and
collaborator changes. This directory is the Git root; its parent contains
separate repositories and software installations.

## Scope and authorization

- Before editing, briefly state the requested outcome, intended changes and
  relevant checks. Keep this proportional to the task.
- Implement only the requested change, using the smallest practical modification
  to the existing structure. Refactoring, renaming, dependency upgrades and new
  constraints require their own authorization when outside the requested scope.
- For example, parallelizing a scan authorizes the changes needed for parallel
  execution; it does not authorize changing the physics model or scan selection.
- Proceed within clear, existing authorization. Ask only when a consequential
  ambiguity or necessary scope expansion remains unresolved. Do not ask again
  for approval already given for the specific change.

## Scientific intent and the comparison model

- The micrOMEGAs model was built on the existing singlet implementation to allow
  comparison in the limit where the TRSM matches the singlet model. Preserving
  that comparison is an intentional scientific requirement.
- Use the collaborators' agreed reference revision, input cards and conventions
  for comparisons. Do not automatically treat the latest implementation as the
  accepted reference. If a required baseline is undocumented, ask the
  collaborators to identify it before making comparison claims.
- Treat model definitions, normalization conventions, effective operators,
  auxiliary fields, numerical prescriptions and scientific defaults as
  protected. Changes, including proposed physics corrections, require explicit
  authorization covering that change.
- This protection includes the authored and generated files in `DM/models/`
  and physics behaviour in `DM/main.c` and `DM/trsm_loop.c`. Read
  `DM/models/README.md` and `DM/models/changes-v2.md` when work touches the model.
  A request to change a wrapper or launcher does not authorize model changes.
- Preserve reference installations, cards and results. Do not silently replace
  a comparison baseline with a corrected or regenerated model.

## Backwards compatibility

- Preserve established commands, configuration formats, defaults, output schemas
  and meanings, and supported checkpoint/resume behaviour unless modification
  is explicitly authorized. A command still running is insufficient if its
  scientific meaning has changed.
- Prefer small additions within the existing code structure. Keep existing
  entry points usable; avoid requiring collaborators to relearn routine tasks.
- If a requested feature requires a compatibility break, explain the concrete
  impact and proposed migration before implementing the break. Keep already
  authorized, independent work moving.

## Findings outside the task

- Finding a suspected bug does not authorize an out-of-scope fix. Report concise
  evidence, confidence, likely impact and a proposed remedy; obtain approval
  before implementing it.
- Raise issues that threaten the requested result promptly. Pause only the work
  whose correctness depends on an unresolved issue.
- Report other relevant findings briefly. Do not turn them into an unsolicited
  audit, implement optional improvements, or repeatedly promote deferred work.

## Focused validation

- Perform checks relevant to the requested change and present compact evidence.
  Reuse existing checks where suitable; do not create a broad verification
  programme for collaborators to complete.
- For parallelization, compare serial and parallel execution on the same small
  set of predefined parameter points, using appropriate numerical tolerances,
  and check for missing or duplicated points.
- For authorized model changes, use the agreed singlet-limit benchmarks and
  checks relevant to the affected physics. Record the model revision, input
  cards, executable identity and comparison tolerances. Distinguish regression
  agreement from independent physics validation.
- Check a small sample before a new production campaign. Preserve running jobs,
  existing outputs and recoverable campaign state; do not interrupt or restart
  them without explicit authorization.
- State what was checked and any material limits. Once the requested work and
  relevant checks are complete, stop rather than expanding the task.

## Documentation and handover

- Update the existing README or walkthrough where appropriate. Keep commands
  directly usable and explanations understandable to a collaborator who has
  not followed the conversation. Avoid unnecessary new documentation files.
- Finish with a short account of what changed, why, how to run it, and what was
  checked. Highlight any authorized change to scientific behaviour or
  compatibility. Support claims with the relevant diff or validation evidence.
