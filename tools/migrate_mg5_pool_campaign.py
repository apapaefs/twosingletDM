#!/usr/bin/env python3
"""Record a validated MG5 scheduling update for a stopped seed campaign.

Every other source, physics setting, executable, checkpoint counter, RNG state,
and saved data prefix must match. Exact old/new code hashes are required.
"""
import argparse
import copy
import importlib.util
import io
import json
from contextlib import redirect_stdout
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import run_trsm_seed_campaign as campaign
from trsm_scan_campaign import CampaignLock, atomic_write_json, configuration_fingerprint, utc_now


def validate_transition(old, new, changes):
    permitted = {"generate_mg5_trsm_xsecs.py", "trsm_constraint_profile.py", "trsm_mg5_pool.py"}
    if not changes or not set(changes).issubset(permitted):
        raise ValueError("Only the MG5 pool implementation/source accounting can change")
    before, after = copy.deepcopy(old), copy.deepcopy(new)
    old_physics, new_physics = before.pop("physics_manifest"), after.pop("physics_manifest")
    old_hashes, new_hashes = old_physics.pop("sources_sha256"), new_physics.pop("sources_sha256")
    old_physics.pop("sha256");new_physics.pop("sha256")
    if old_physics != new_physics:
        raise ValueError("Physics manifest changed beyond source accounting")
    actual = {k: {"before": old_hashes.get(k), "after": new_hashes.get(k)}
              for k in set(old_hashes) | set(new_hashes) if old_hashes.get(k) != new_hashes.get(k)}
    if actual != changes:
        raise ValueError(f"Source changes differ from the reviewed transition: {actual}")
    pool = after["mg5_runtime"].pop("pool", None)
    if not pool or before != after:
        raise ValueError("Runtime changed beyond the validated MG5 pool")
    return pool


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--source-transition", type=Path, required=True)
    parser.add_argument("--validation-report", type=Path, required=True)
    args = parser.parse_args()
    directory = args.campaign_dir.resolve()
    state_path = directory / "campaign_state.json"
    changes = json.loads(args.source_transition.read_text())["source_changes"]
    validation = json.loads(args.validation_report.read_text())
    if validation.get("passed") is not True:
        raise ValueError("A passing native serial/pooled validation report is required")
    with CampaignLock(state_path, sys.argv):
        state = json.loads(state_path.read_text())
        if configuration_fingerprint(state["configuration"]) != state["configuration_fingerprint"]:
            raise ValueError("Saved campaign configuration fingerprint does not match")
        campaign.assert_no_live_workers(state)
        if state["status"] == "running":
            raise ValueError("Stop the campaign cleanly before migrating it")
        cargs = campaign.parse_args(["--campaign-dir", str(directory), "--resume"])
        campaign.hydrate_configuration(cargs, state)
        receipt = campaign.run_preflight(cargs)
        pool = validate_transition(state["runtime_receipt"], receipt, changes)
        if validation.get("pool_root") != pool["root"]:
            raise ValueError("Native validation report belongs to another pool")
        native_rows = validation.get("serial", []) + validation.get("parallel", [])
        expected_driver = changes["generate_mg5_trsm_xsecs.py"]["after"]
        if not native_rows or any(row.get("driver_sha256") != expected_driver for row in native_rows):
            raise ValueError("Native validation used a different MG5 driver")
        first = state["seeds"][0]
        generator = Path(state["configuration"]["generator_script"])
        previous_argv = sys.argv
        command = first["command"]
        index = command.index(str(generator))
        sys.argv = [str(generator), *command[index+1:]]
        try:
            spec = importlib.util.spec_from_file_location("_pool_migration_generator", generator)
            g = importlib.util.module_from_spec(spec)
            with redirect_stdout(io.StringIO()):spec.loader.exec_module(g)
        finally:
            sys.argv = previous_argv
        updates = []
        for record in state["seeds"]:
            checkpoint = campaign.seed_checkpoint(record)
            if checkpoint is None:
                raise ValueError(f"Missing checkpoint for seed {record['seed']}")
            campaign.validate_prefix(record["point_output"], checkpoint["outputs"]["main"])
            point = Path(record["point_output"])
            metadata_path = point.with_suffix(".metadata.json")
            metadata = json.loads(metadata_path.read_text())
            options = g.parse_args(["--resume-from", str(point), "--nrandom", str(cargs.nrandom)])
            if options.planck_cmb:
                options._cmb_driver = g.require_cmb_capability(Path(g.micromegas_configuration(options)["executable"]))
            new_config = g.immutable_scan_configuration(options)
            old_config = copy.deepcopy(new_config)
            old_config["physics_manifest"] = metadata["physics_manifest"]
            old_fingerprint = configuration_fingerprint(old_config)
            if old_fingerprint != checkpoint["configuration_fingerprint"] or old_fingerprint != metadata["configuration_fingerprint"]:
                raise ValueError(f"Saved options/physics do not match seed {record['seed']}")
            if metadata["physics_manifest"] != state["runtime_receipt"]["physics_manifest"]:
                raise ValueError(f"Seed {record['seed']} has different physics provenance")
            if new_config["physics_manifest"] != receipt["physics_manifest"]:
                raise ValueError("Preflight and generator manifests differ")
            new_fingerprint = configuration_fingerprint(new_config)
            checkpoint["configuration_fingerprint"] = new_fingerprint
            metadata["configuration_fingerprint"] = new_fingerprint
            metadata["physics_manifest"] = receipt["physics_manifest"]
            metadata.setdefault("operational_updates", []).append({"updated_utc": utc_now(),
                "kind": "independent_mg5_pool", "pool": pool,
                "old_configuration_fingerprint": old_fingerprint,
                "new_configuration_fingerprint": new_fingerprint})
            updates.extend([(point.with_suffix(".checkpoint.json"), checkpoint), (metadata_path, metadata)])
        backup = directory / "operational-updates" / "mg5-pool"
        backup.mkdir(parents=True, exist_ok=False)
        report = {"schema": "trsm_mg5_pool_transition_v1", "created_utc": utc_now(),
                  "source_changes": changes, "validation": validation,
                  "old_runtime_receipt": state["runtime_receipt"], "new_runtime_receipt": receipt,
                  "seed_count": len(state["seeds"]), "pool": pool,
                  "data_files_changed": False}
        atomic_write_json(backup / "transition.json", report)
        # Back up every file before modifying any of them.
        for path in [state_path, *(path for path, _ in updates)]:
            saved = backup / "before" / path.relative_to(directory)
            saved.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, saved)
        for path, payload in updates:
            atomic_write_json(path, payload)
        state["runtime_receipt"] = receipt
        state["runtime_fingerprint"] = configuration_fingerprint(receipt)
        state.setdefault("operational_updates", []).append(str(backup / "transition.json"))
        atomic_write_json(state_path, state)
        print(json.dumps({"migrated_seeds": len(state["seeds"]), "backup": str(backup), "pool": pool}, indent=2))


if __name__ == "__main__":
    main()
