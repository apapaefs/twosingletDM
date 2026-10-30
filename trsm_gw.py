"""BSMPT 3.2.1 CalcGW adapter; keep quantities paired by transition index."""

import json
import math
import os
from pathlib import Path
import re
import subprocess

from test_trsm_ewpt import (
    DEFAULT_EXECUTABLE, TRSMEWPTPoint, parse_calctemps_output, write_trsm_input,
)
from trsm_inputs import json_safe

DEFAULT_GW_EXECUTABLE = DEFAULT_EXECUTABLE.with_name("CalcGW")
GW_FIELDS = {
    "T_star_GeV": "T_star", "T_reh_GeV": "T_reh", "v_wall": "v_wall",
    "alpha": "alpha_PT", "beta_over_H": "beta/H",
    "kappa_col": "kappa_col", "kappa_sw": "kappa_sw", "eps_turb": "eps_turb",
    "cs_false": "cs_f", "cs_true": "cs_t",
    "f_col_Hz": "fb_col", "h2omega_col": "h2Omegab_col",
    "f1_sw_Hz": "f_1_sw", "f2_sw_Hz": "f_2_sw", "h2omega2_sw": "h2Omega_2_sw",
    "f1_turb_Hz": "f_1_turb", "f2_turb_Hz": "f_2_turb", "h2omega2_turb": "h2Omega_2_turb",
    "snr_lisa_col": "SNR(LISA-3yrs)_col", "snr_lisa_sw": "SNR(LISA-3yrs)_sw",
    "snr_lisa_turb": "SNR(LISA-3yrs)_turb", "snr_lisa": "SNR(LISA-3yrs)",
}
GW_COLUMNS = ("gw_status", "gw_error", "gw_native_summary", "gw_transitions",
              "gw_transition_count", "gw_successful_transition_count", "gw_output_path")


def finite(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def parse_gw_row(row, temperature_kind="perc"):
    """Keep failed native transitions too, but never promote them to predictions."""
    required = ("status_nlo_stability", "status_tracing", "status_coex_pairs")
    if any(key not in row for key in required):
        raise ValueError("CalcGW output lacks native assessment statuses")
    summary = {key: row.get(key) for key in (*required, "status_ewsr", "transition_history")}
    base_success = all(row[key] == "success" for key in required)
    indices = sorted(int(m[1]) for key in row if (m := re.fullmatch(r"status_gw_(\d+)", key)))
    transitions = []
    for index in indices:
        transition = {"transition_index": index, "temperature_kind": temperature_kind,
                      "status": row[f"status_gw_{index}"],
                      "temperature_status": row.get(f"status_{temperature_kind}_{index}")}
        transition.update({name: finite(row.get(f"{native}_{index}"))
                           for name, native in GW_FIELDS.items()})
        transition["valid"] = (
            base_success and transition["status"] == "success"
            and transition["temperature_status"] == "success"
            and all(transition[key] is not None and transition[key] > 0
                    for key in ("T_star_GeV", "alpha", "beta_over_H"))
            and transition["beta_over_H"] >= 1
        )
        transitions.append(transition)
    if row["status_coex_pairs"] == "success" and not indices:
        raise ValueError("CalcGW output has coexistence pairs but no GW fields")
    status = "success" if any(t["valid"] for t in transitions) else "no_successful_gw"
    if row["status_nlo_stability"] != "success" or row["status_tracing"] != "success" \
            or row["status_coex_pairs"] not in ("success", "no_coex_pair"):
        status = "native_failure"
    return json_safe({"gw_status": status, "gw_error": "", "gw_native_summary": summary,
            "gw_transitions": transitions, "gw_transition_count": len(transitions),
            "gw_successful_transition_count": sum(t["valid"] for t in transitions)})


def run_calcgw(row, index, directory, options):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    input_path, output_path = directory / "TRSM_Input.tsv", directory / "CalcGW.tsv"
    result = {"gw_status": "invalid_input", "gw_error": "", "gw_native_summary": {},
              "gw_transitions": [], "gw_transition_count": 0,
              "gw_successful_transition_count": 0, "gw_output_path": str(output_path)}
    try:
        keys = {"m2": "M2", "m3": "M3", "vs": "vs", "a12": "a12",
                "lx": "lX", "lphix": "lPhiX", "lsx": "lSX"}
        values = {name: finite(row.get(column)) for name, column in keys.items()}
        if any(value is None for value in values.values()):
            raise ValueError("missing/nonfinite TRSM input")
        if finite(row.get("vx")) != 0 or any(values[key] <= 0 for key in ("m2", "m3")) or values["vs"] == 0:
            raise ValueError("CalcGW adapter requires vx=0, nonzero vs, and positive masses")
        write_trsm_input(input_path, TRSMEWPTPoint(index=index, **values))
        command = [str(options["executable"]), "--model=trsm", f"--input={input_path}",
                   f"--output={output_path}", "--firstline=2", "--lastline=2",
                   "--multistepmode=default", "--usemultithreading=false",
                   *(f"--{key}={options[key]}" for key in
                     ("thigh", "trans_temp", "vwall", "epsturb"))]
        # Each invocation is isolated and single threaded; the launcher bounds jobs.
        (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
        output_path.unlink(missing_ok=True)
        environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        with (directory / "CalcGW.log").open("w") as log:
            result["gw_status"] = "execution_failed"
            completed = subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT,
                                       timeout=options["timeout"], check=False, env=environment)
        if completed.returncode:
            raise ValueError(f"CalcGW exited with code {completed.returncode}; see CalcGW.log")
        result["gw_status"] = "invalid_output"
        rows = parse_calctemps_output(output_path)
        if len(rows) != 1 or finite(rows[0].get("index")) != index:
            raise ValueError("CalcGW must return exactly the requested point")
        result.update(parse_gw_row(rows[0], options["trans_temp"]))
    except subprocess.TimeoutExpired:
        result.update(gw_status="timeout", gw_error=f"CalcGW exceeded {options['timeout']:g} seconds")
    except (OSError, ValueError) as error:
        result["gw_error"] = str(error)
    return result
