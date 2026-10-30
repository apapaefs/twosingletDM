#!/usr/bin/env python3
"""Inspect a v2 TRSM TSV without installing the simulation software."""

import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import sys


OVERVIEW = [
    "seed", "point_index", "M2", "M3", "vs", "a12", "k1", "k2",
    "lX", "lPhiX", "lSX", "K133", "K233", "w1", "w2",
    "h1_h3h3_br", "h2_h3h3_br", "h2_h1h1_br", "dm_omega",
    "dm_dir_det", "dm_dir_det_limit", "dm_indirect_ratio", "dm_cmb_ratio",
    "evo", "thc", "experimental_subset", "flavour", "dm",
    "vacuum_tree_global", "rg_bfb", "rg_unitarity", "theory_strict_subset",
    "ewpt_status", "ewpt_ew_entry_jump_over_T", "ewpt_ew_entry_temperature_GeV",
    "ewpt_baryo_candidate", "ewpt_gw_candidate",
    "mg5_xsec_gg_heta0_pb", "mg5_xsec_pp_eta0Z_total_pb", "mg5_xsec_hh_pb",
    "mono_higgs_xsec_pb", "mono_z_total_xsec_pb", "gw_status",
    "gw_successful_transition_count",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Original TSV or gzip-compressed TSV")
    parser.add_argument("--out", type=Path, default=Path("quicklook"))
    args = parser.parse_args()
    # The table contains JSON cells, which can exceed the csv default limit.
    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    opener = gzip.open if args.input.suffix == ".gz" else open
    counts = Counter(rows=0, evo_and_thc=0, experimental=0,
                     non_dm_viable=0, fully_viable=0,
                     calcgw_successful_points=0, fully_viable_calcgw_points=0)
    args.out.mkdir(parents=True, exist_ok=True)
    paths = [args.out / "fully_viable_points.tsv.gz",
             args.out / "fully_viable_overview.csv", args.out / "summary.json"]
    if any(path.exists() for path in paths):
        parser.error("Output files already exist; choose a new --out directory.")
    with opener(args.input, "rt", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        fields = reader.fieldnames or []
        required = {"seed", "point_index", "M2", "M3", "constraint_version",
                    "thc", "experimental_subset", "flavour", "dm", "evo"}
        if required - set(fields):
            parser.error(f"Missing columns: {sorted(required - set(fields))}")
        overview = [name for name in OVERVIEW if name in fields]
        with gzip.open(paths[0], "wt", newline="") as full_stream, \
                paths[1].open("w", newline="") as overview_stream:
            full_writer = csv.DictWriter(full_stream, fieldnames=fields, delimiter="\t")
            simple_writer = csv.DictWriter(overview_stream, fieldnames=overview,
                                          extrasaction="ignore")
            full_writer.writeheader()
            simple_writer.writeheader()
            for row in reader:
                if row["constraint_version"] != "trsm_constraints_v2":
                    raise ValueError("This reader requires trsm_constraints_v2 throughout.")
                # Only an explicit True supplies a pass; nan means unassessed.
                experimental = row["experimental_subset"] == "True"
                non_dm = row["thc"] == "True" and experimental and row["flavour"] == "True"
                full = non_dm and row["dm"] == "True"
                gw = row.get("gw_status") == "success"
                counts["rows"] += 1
                counts["evo_and_thc"] += row["evo"] == "True" and row["thc"] == "True"
                counts["experimental"] += experimental
                counts["non_dm_viable"] += non_dm
                counts["fully_viable"] += full
                counts["calcgw_successful_points"] += gw
                counts["fully_viable_calcgw_points"] += full and gw
                if full:
                    full_writer.writerow(row)
                    simple_writer.writerow(row)
    summary = {"input": str(args.input.resolve()), "counts": dict(counts),
               "fully_viable_selection": "thc & experimental_subset & flavour & dm",
               "required_value": "True", "point_key": ["seed", "point_index"]}
    paths[2].write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"\nWritten to {args.out.resolve()}:")
    for path in paths:
        print(f"  {path.name}")


if __name__ == "__main__":
    main()
