#!/usr/bin/env python3
"""Export the installed HiggsTools analysis catalogue without evaluating points.

The scope labels identify unavailable model topologies, not per-point acceptance.
Saved selected-limit counts can be included without rerunning HiggsBounds.
"""

import argparse
from collections import Counter
import importlib.metadata
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from trsm_constraint_profile import dataset_digest, git_revision
from trsm_paths import higgs_dataset_path


def scope(data, relative_path):
    """Only rule out channels known to vanish in the current vx=0 adapter."""
    if relative_path.parts[0] != "neutral":
        return "C", "No charged scalar in this model."
    if '"CPValue": "odd"' in json.dumps(data.get("constraints", {})):
        return "P", "Requires a CP-odd scalar; both visible scalars are CP even."
    if data["limitClass"] == "ChainDecayLimit":
        return "Z", "No scalar-to-Z-scalar decay in this adapter."
    if data["limitClass"] == "PairProductionLimit" and data["collider"] == "LEP":
        return "Z", "No LEP scalar-pair production input."
    processes = data["process"]
    processes = processes if isinstance(processes, list) else [processes]
    channels = [c for p in processes for c in p.get("channels", [])]
    decays = {"uu", "dd", "ee", "bb", "tautau", "mumu", "cc", "ss", "tt", "gg", "gamgam", "Zgam", "WW", "ZZ", "inv", "none"}
    if channels and all(p in {"brtHu", "brtHc"} or d not in decays for p, d in channels):
        return "Z", "Requires a production or decay input that is zero in this adapter."
    return "A", "Potentially applicable; the provider checks the point's masses, widths and analysis conditions."


def code(value):
    return r"\code{" + str(value) + "}"


def compact_process(data, native):
    """Group repeated Cartesian channel lists without changing their meaning."""
    if len(native) <= 300 or not isinstance(data["process"], dict):
        return native
    channels = data["process"].get("channels", [])
    if not channels:
        return native
    by_production = {}
    for production, decay in channels:
        by_production.setdefault(production, set()).add(decay)
    grouped = {}
    for production, decays in by_production.items():
        grouped.setdefault(tuple(sorted(decays)), []).append(production)
    return data["collider"] + " " + " + ".join(
        "[" + ",".join(productions) + "]>[" + ",".join(decays) + "]"
        for decays, productions in grouped.items())


def reference(value, filename, dataset, commit):
    if re.fullmatch(r"(?:\d{4}\.\d{4,5}|hep-ex/\d{7})", value):
        return r"\href{https://arxiv.org/abs/" + value + "}{" + code(value) + "}"
    return (r"\href{https://gitlab.com/higgsbounds/" + dataset + "/-/blob/" + commit
            + "/" + filename + "}{" + code(value or "dataset entry") + "}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "docs/higgstools-analyses.tex")
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--selected-counts", type=Path)
    parser.add_argument("--selected-only", action="store_true",
                        help="Show only HiggsBounds rows with a positive saved selection count.")
    args = parser.parse_args(argv)
    if args.selected_only and args.selected_counts is None:
        parser.error("--selected-only requires --selected-counts")
    import Higgs.bounds as HB
    import Higgs.signals as HS

    selected = json.loads(args.selected_counts.read_text()) if args.selected_counts else {}
    counts = selected.get("selected_point_counts", {})
    inventory = {"higgstools_version": importlib.metadata.version("HiggsTools"),
                 "datasets": {}, "bounds": [], "signals": [], "selected_counts_source": selected}
    for name in ("hbdataset", "hsdataset"):
        path = higgs_dataset_path(name)
        inventory["datasets"][name] = {"commit": git_revision(path), "tree_sha256": dataset_digest(path)}
    bounds = HB.Bounds(str(higgs_dataset_path("hbdataset")))
    for limit in bounds.limits():
        path = Path(limit.loadedFrom())
        data = json.loads(path.read_text())
        relative = path.resolve().relative_to(higgs_dataset_path("hbdataset").resolve())
        label, reason = scope(data, relative)
        inventory["bounds"].append({
            "id": limit.id(), "reference": limit.reference(), "cite_key": limit.citeKey(),
            "process": limit.processDesc(), "display_process": compact_process(data, limit.processDesc()),
            "domain": limit.extentDesc(),
            "experiment": data["experiment"], "collider": data["collider"],
            "luminosity_fb": limit.luminosity(), "file": str(relative),
            "scope": label, "scope_reason": reason, "constraints": data.get("constraints", {}),
            "selected_points": counts.get(str(limit.id()), 0) if selected else None,
        })
    signals = HS.Signals(str(higgs_dataset_path("hsdataset")))
    for measurement in signals.measurements():
        path = Path(measurement.loadedFrom())
        data = json.loads(path.read_text())
        relative = path.resolve().relative_to(higgs_dataset_path("hsdataset").resolve())
        inventory["signals"].append({
            "id": measurement.id(), "reference": measurement.reference(),
            "cite_key": measurement.citeKey(), "experiment": data["experiment"],
            "collider": data["collider"], "luminosity_fb": measurement.luminosity(),
            "file": str(relative), "submeasurements": len(measurement.subMeasurements()),
            "processes": sorted({m.processDesc() for m in measurement.subMeasurements().values()}),
        })
    for key in ("bounds", "signals"):
        assert len({x["id"] for x in inventory[key]}) == len(inventory[key])
        inventory[key].sort(key=lambda x: (x.get("scope", ""), x["file"]))
    unknown_ids = set(counts) - {str(x["id"]) for x in inventory["bounds"]}
    if unknown_ids:
        raise ValueError(f"Saved selected IDs absent from dataset: {sorted(unknown_ids)}")
    wrong_scope = [x["id"] for x in inventory["bounds"] if x["selected_points"] and x["scope"] != "A"]
    if wrong_scope:
        raise ValueError(f"Scope disagrees with saved selections: {wrong_scope}")
    args.inventory.parent.mkdir(parents=True, exist_ok=True)
    args.inventory.write_text(json.dumps(inventory, indent=2) + "\n")

    scope_counts = Counter(x["scope"] for x in inventory["bounds"])
    table_bounds = [x for x in inventory["bounds"]
                    if not args.selected_only or x["selected_points"] > 0]
    lines = [r"% Generated by tools/export_higgstools_tables.py; do not edit rows manually.",
             r"\section{HiggsTools analyses in the installed datasets}",
             r"\label{sec:catalogue}",
             "The installed datasets contain "
             + f"{len(inventory['bounds'])} HiggsBounds limits and {len(inventory['signals'])} HiggsSignals measurement files. "
             + (f"The HiggsBounds table lists the {len(table_bounds)} IDs with a nonzero selection count in the reference campaign; "
                "rows with a zero count are omitted. All HiggsSignals measurement files are listed. "
                if args.selected_only else "The tables below list all of these entries. ")
             + "Each row is one provider ID; several rows may refer to the same experimental paper. "
             + "These are the datasets used by the reference campaign, rather than a claim to include later experimental results.",
             r"\begin{itemize}"]
    for name, record in inventory["datasets"].items():
        lines.append(r"\item " + code(name) + " revision " + code(record["commit"]) + r".\newline SHA256: " + code(record["tree_sha256"]) + ".")
    lines += [r"\item Scope A: the input topology can be nonzero in this adapter. "
              "Mass, width, kinematic and model-likeness conditions still determine applicability at each point; "
              "an A label does not assert sensitivity or exclusion throughout either scan rectangle."]
    if not args.selected_only:
        lines += [r"\item Scope C: charged-scalar analyses; P: CP-odd-only analyses; "
                  "Z: a required production or decay input is zero. These entries are listed separately for completeness."]
    lines += [r"\item The domain is HiggsTools' own grid description. Masses are in GeV; "
              "width coordinates and other qualifications follow the linked implementation. "
              "Long decimal endpoints are displayed to six significant figures; the exported JSON catalogue keeps the original strings. "
              "The table does not extrapolate a limit outside its applicability region.",
              r"\item Process notation follows HiggsTools: H denotes a neutral scalar; "
              "ggH and vbfH denote gluon fusion and vector-boson fusion; HW, HZ, Htt and bbH denote associated production; "
              "eeHZ denotes LEP Higgsstrahlung; inv denotes invisible decay; gam denotes a photon. "
              "The symbol > separates production and decay, and comma-separated entries form an analysis combination.",
              r"\item X1, X2 and X3 and the corresponding M1, M2 and M3 in the provider descriptions "
              "label the parent and daughter roles of that search. They are not fixed identifications with "
              "the model's $h_1,h_2,h_3$ or their masses. HiggsBounds tests allowed particle assignments.",
              r"\end{itemize}"]
    if selected:
        lines += ["The final column counts distinct reference-campaign points for which the analysis was "
                  "selected for at least one scalar by expected sensitivity, using all "
                  + f"{selected['rows']:,} saved rows. A zero count does not mean a search is inapplicable; "
                  "only selected limits, not every tested limit, are recorded exhaustively. "
                  + f"There are {sum(v > 0 for v in counts.values())} different selected IDs."]
    for group, title in (("A", "HiggsBounds searches selected in the reference campaign"
                          if args.selected_only else "Potentially applicable HiggsBounds searches"),
                         ("CPZ", "Other loaded HiggsBounds searches")):
        records = [x for x in table_bounds if x["scope"] in group]
        if not records:
            continue
        lines += [r"\begin{landscape}", r"\subsection{" + title + "}", r"\footnotesize",
                  r"\setlength{\tabcolsep}{3pt}",
                  r"\begin{longtable}{@{}p{21mm}p{75mm}p{48mm}p{27mm}p{29mm}r@{}}",
                  r"\toprule ID / scope & Search process & Provider domain & Experiment / $\mathcal{L}$ [fb$^{-1}$] & Reference & Selected\\\midrule\endhead"]
        for x in records:
            ref = reference(x["reference"], x["file"], "hbdataset", inventory["datasets"]["hbdataset"]["commit"])
            display_domain = re.sub(r"\d+\.\d+(?:e[+-]?\d+)?", lambda m: f"{float(m[0]):.6g}", x["domain"])
            domain = r"\newline ".join(r"\mbox{" + code(part.strip()) + "}"
                                       for part in re.split(r"(?<=\)),", display_domain))
            cells = [code(x["id"]) + r"\newline " + x["scope"], code(x["display_process"]), domain,
                     x["experiment"] + r"\newline " + f"{x['luminosity_fb']:g}", ref,
                     f"{x['selected_points']:,}" if selected else "---"]
            lines.append(" & ".join(cells) + r"\\")
        lines += [r"\bottomrule\end{longtable}\end{landscape}\normalsize"]
    lines += [r"\begin{landscape}", r"\subsection{HiggsSignals measurements}",
              "All files below enter the HiggsSignals likelihood. The submeasurement count identifies "
              "the number of channels or bins in each file; correlations are supplied by the dataset.",
              r"\small",
              r"\begin{longtable}{@{}p{25mm}p{100mm}p{35mm}p{32mm}r@{}}",
              r"\toprule ID & Measurement file / channel & Experiment / $\mathcal{L}$ [fb$^{-1}$] & Reference & Bins\\\midrule\endhead"]
    for x in inventory["signals"]:
        ref = reference(x["reference"], x["file"], "hsdataset", inventory["datasets"]["hsdataset"]["commit"])
        lines.append(" & ".join([code(x["id"]), code(x["file"]),
                                 x["experiment"] + " " + x["collider"] + r"\newline " + f"{x['luminosity_fb']:g}",
                                 ref, str(x["submeasurements"])]) + r"\\")
    lines += [r"\bottomrule\end{longtable}\end{landscape}\normalsize"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n")
    print(json.dumps({"HB": len(inventory["bounds"]), "HS": len(inventory["signals"]),
                      "HB_table_rows": len(table_bounds),
                      "scope_counts": dict(scope_counts), "output": str(args.output)}))


if __name__ == "__main__":
    main()
