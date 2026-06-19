import argparse
import json
import math
import os
import random
from datetime import date
from pathlib import Path


TRSM_POINT_ARGS = ("m2", "m3", "vs", "a12", "lx", "lphix", "lsx")
EWPT_STRENGTH_PRIORITY = ("nucl", "perc", "compl", "crit")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Generate TRSM points. By default this runs the existing random "
            "vx=0 scan; providing all explicit point parameters skips the "
            "random scan and evaluates that one point."
        )
    )
    parser.add_argument(
        "seed",
        nargs="?",
        type=int,
        default=1,
        help="Random seed for the default random scan, or bookkeeping seed for a single point.",
    )
    parser.add_argument(
        "--m1",
        type=float,
        default=125.09,
        help="Accepted for consistency with test_trsm_ewpt.py; generate_lams currently fixes M1=125.09.",
    )
    parser.add_argument("--m2", type=float)
    parser.add_argument("--m3", type=float)
    parser.add_argument("--vs", type=float)
    parser.add_argument("--a12", type=float)
    parser.add_argument("--lx", type=float)
    parser.add_argument("--lphix", type=float)
    parser.add_argument("--lsx", type=float)
    resonant_group = parser.add_mutually_exclusive_group()
    resonant_group.add_argument(
        "--resonantDM1",
        action="store_true",
        help="Set m3 = 2*m1 for every evaluated point.",
    )
    resonant_group.add_argument(
        "--resonantDM2",
        action="store_true",
        help="Set m3 = 2*m2 for every evaluated point.",
    )
    parser.add_argument(
        "--nrandom",
        type=int,
        default=100,
        help="Number of points for the default random scan.",
    )
    parser.add_argument(
        "--scan-k133-k233",
        action="store_true",
        help=(
            "In the random vx=0 scan, sample K133 and K233 instead of "
            "lambda_phiX and lambda_SX, then derive lphix/lsx with the "
            "current coupling convention."
        ),
    )
    parser.add_argument(
        "--write-dm-failed",
        action="store_true",
        help=(
            "Also write points that pass all non-DM vx=0 checks but fail the "
            "dark-matter check to trsm_points_<tag>_dm_failed.dat."
        ),
    )
    parser.add_argument(
        "--write-all-points",
        action="store_true",
        help=(
            "Write every evaluated vx=0 point to the main output file, "
            "regardless of constraints, and run EWPT for every point."
        ),
    )
    parser.add_argument(
        "--print-info",
        action="store_true",
        help=(
            "Print the full vx=0 point table, constraint summary, and dark-matter "
            "summary for each fully evaluated scan point. This does not change "
            "which points are written."
        ),
    )
    parser.add_argument(
        "--run-ewpt",
        action="store_true",
        help="Run test_trsm_ewpt.py/BSMPT only after a point passes all viability checks.",
    )
    parser.add_argument(
        "--run-ewpt-on-dm-failed",
        action="store_true",
        help=(
            "Also run test_trsm_ewpt.py/BSMPT for vx=0 points that pass all "
            "non-DM checks but fail the dark-matter check. These points are "
            "written to the _dm_failed sidecar."
        ),
    )
    parser.add_argument(
        "--higgstools-details",
        action="store_true",
        help=(
            "Print the most constraining HiggsBounds limits and largest "
            "HiggsSignals chi-squared contributors for each evaluated point."
        ),
    )
    parser.add_argument(
        "--higgstools-top",
        type=int,
        default=5,
        help="Number of HiggsTools analyses/contributors to print with --higgstools-details.",
    )
    parser.add_argument(
        "--no-save-higgstools-details",
        dest="save_higgstools_details",
        action="store_false",
        default=True,
        help="Do not save compact HiggsBounds/HiggsSignals limiting-analysis summaries to the scan output.",
    )
    parser.add_argument(
        "--ewpt-executable",
        type=Path,
        help="CalcTemps executable passed to the EWPT runner.",
    )
    parser.add_argument(
        "--ewpt-minima-executable",
        type=Path,
        help="MinimaTracer executable passed to the EWPT runner.",
    )
    parser.add_argument(
        "--ewpt-thigh",
        type=float,
        default=300.0,
        help="High-temperature bound passed to MinimaTracer and CalcTemps.",
    )
    parser.add_argument(
        "--ewpt-multistepmode",
        default="default",
        help="BSMPT multistepmode passed to MinimaTracer and CalcTemps.",
    )
    parser.add_argument(
        "--ewpt-workdir",
        type=Path,
        help="Base directory for EWPT outputs; each viable point gets a point_XXXXXX subdirectory.",
    )
    parser.add_argument(
        "--ewpt-plot-phases",
        action="store_true",
        help="Write MinimaTracer phase plots for viable points.",
    )
    parser.add_argument(
        "--ewpt-plot-output",
        type=Path,
        help="Phase plot basename inside each viable point EWPT workdir.",
    )
    parser.add_argument(
        "--ewpt-plot-format",
        choices=["png", "pdf", "both"],
        default="both",
        help="Phase plot format for --ewpt-plot-phases.",
    )
    parser.add_argument(
        "--ewpt-require-eq418",
        action="store_true",
        help="Before running EWPT, require the Eq. 4.18 quartic-coupling positivity check.",
    )
    parser.add_argument("--ewpt-sym-threshold", type=float, default=1.0)
    parser.add_argument("--ewpt-w1-threshold", type=float, default=5.0)
    parser.add_argument("--ewpt-wx-threshold", type=float, default=1.0)
    parser.add_argument("--ewpt-ws-threshold", type=float, default=1.0)
    args = parser.parse_args(argv)

    provided = [name for name in TRSM_POINT_ARGS if getattr(args, name) is not None]
    resonance_enabled = args.resonantDM1 or args.resonantDM2
    required_point_args = [
        name
        for name in TRSM_POINT_ARGS
        if name != "m3" or not resonance_enabled
    ]
    if provided and any(getattr(args, name) is None for name in required_point_args):
        missing = ", ".join(
            f"--{name}" for name in required_point_args if getattr(args, name) is None
        )
        parser.error(f"explicit point mode requires all point parameters; missing {missing}")
    if provided and not math.isclose(args.m1, 125.09, rel_tol=0.0, abs_tol=1e-9):
        parser.error("--m1 is accepted for CLI compatibility, but this generator currently fixes M1=125.09")
    if args.higgstools_top < 1:
        parser.error("--higgstools-top must be at least 1")
    return args


cli_args = parse_args()
ini_seed = cli_args.seed

from generate_trsm_info import * # TRSM info generator (branching ratios, mixing matrices, etc.)
from test_trsm_evolution import * # RGE evolution
#from twosinglet_sigmahhh_read_pickle_nn import * # Triple Higgs Cross Section
from test_trsm_higgstools import * # HiggsTools setup 
from trsm_kstoalphas import * # Convert from k1, k2, k3 to a12, a13, a23
from generate_mg5_trsm_xsecs import * # call MG5 to get the cross section for a specific proces. Make sure that the process has been generated (and check run card for energy/cuts etc!)
from mg5_process_runner import run_mg5_processes # run selected MG5 processes and collect cross sections
from scan_output import write_valid_point as write_valid_point_file
from test_trsm_theory_constraints import * # unitarity/boundedness from below
from test_trsm_DM import print_dm_info, test_dm # micrOMEGAs dark matter check for the vx=0 branch
from prettytable import PrettyTable
from singlet_EWPO import * # Electroweak Precision Observables

#############
# OPTIONS
##############

# print debug?
debug = False

# run MG5 on points that pass constraints?
RunMG5 = False

# MG5 processes to run when RunMG5 is True.
# Available process names are defined by ProcLocation in generate_mg5_trsm_xsecs.py.
MG5ProcessesToRun = [] # e.g. ['hh', 'hhh']

# for random scan within ranges, how many points to run
nrandom=cli_args.nrandom

# Energy for xsec calculation
Energy = 13.6

# HH SM XSEC at Energy [pb] (nn23lo1 pdf in MadGraph)
xsec_sm = {}
xsec_sm[13] = 0.01452
xsec_sm[13.6] = 0.01617

# set the Higgs mass:
mhiggs = 125.

# TAG for RUN output
RunTag = str(Energy) + '-' + str(date.today()).replace('-','') + '-' + str(ini_seed) + '-' + str(RunMG5)

# Directory for output:
OutputDir = 'output/'

# reset the output before next run?
ResetOutput = True

# Print additional TRSM point info?
PRINTINFO = False

SM_VEV_FOR_K_SCAN = 246.

###########################################################
# some functions
###########################################################

def lambdas_to_k133_k233(lPhiX, lSX, vs, a12, v=SM_VEV_FOR_K_SCAN):
    c = math.cos(a12)
    s = math.sin(a12)
    k133 = 0.5 * (lPhiX * v * c - lSX * vs * s)
    k233 = 0.5 * (lPhiX * v * s + lSX * vs * c)
    return k133, k233


def k133_k233_to_lambdas(K133, K233, vs, a12, v=SM_VEV_FOR_K_SCAN):
    c = math.cos(a12)
    s = math.sin(a12)
    lPhiX = 2.0 * (c * K133 + s * K233) / v
    lSX = 2.0 * (-s * K133 + c * K233) / vs
    return lPhiX, lSX

# print the parameter point info info:
def print_info(vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3):
    
    tbl = PrettyTable(["var", "value"])
    tbl.add_row(['vs', vs])
    tbl.add_row(['vx', vx])
    tbl.add_row(['a12', a12])
    tbl.add_row(['a13', a13])
    tbl.add_row(['a23', a23])
    tbl.add_row(['c1', math.cos(a12)])
    tbl.add_row(['c2', math.cos(a13)])
    tbl.add_row(['c3', math.cos(a23)])
    tbl.add_row(['s1', math.sin(a12)])
    tbl.add_row(['s2', math.sin(a13)])
    tbl.add_row(['s3', math.sin(a23)])
    tbl.add_row(['M2', M2])
    tbl.add_row(['w2', w2])
    tbl.add_row(['M3', M3])
    tbl.add_row(['w3', w3])
    tbl.add_row(['k1', k1])
    tbl.add_row(['k2', k2])
    tbl.add_row(['k3', k3])
    tbl.add_row(['K111', K111])
    tbl.add_row(['K112', K112])
    tbl.add_row(['K113', K113])
    tbl.add_row(['K123', K123])
    tbl.add_row(['K122', K122])
    tbl.add_row(['K1111', K1111])
    tbl.add_row(['K1112', K1112])
    tbl.add_row(['K1113', K1113])
    tbl.add_row(['K133', K133])
    print(tbl)
    #print('\n')

# print the parameter point info for vx=0:
def print_info_vxzero(vs, vx, M2, M3, a12, a13, a23, lX, lPhiX, lSX, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3):
    
    tbl = PrettyTable(["var", "value"])
    tbl.add_row(['vs', vs])
    tbl.add_row(['vx', vx])
    tbl.add_row(['a12', a12])
    tbl.add_row(['a13', a13])
    tbl.add_row(['a23', a23])
    tbl.add_row(['c1', math.cos(a12)])
    tbl.add_row(['c2', math.cos(a13)])
    tbl.add_row(['c3', math.cos(a23)])
    tbl.add_row(['s1', math.sin(a12)])
    tbl.add_row(['s2', math.sin(a13)])
    tbl.add_row(['s3', math.sin(a23)])
    tbl.add_row(['lX', lX])
    tbl.add_row(['lPhiX', lPhiX])
    tbl.add_row(['lSX', lSX])
    tbl.add_row(['M2', M2])
    tbl.add_row(['w2', w2])
    tbl.add_row(['M3', M3])
    tbl.add_row(['w3', w3])
    tbl.add_row(['k1', k1])
    tbl.add_row(['k2', k2])
    tbl.add_row(['k3', k3])
    tbl.add_row(['K111', K111])
    tbl.add_row(['K112', K112])
    tbl.add_row(['K113', K113])
    tbl.add_row(['K123', K123])
    tbl.add_row(['K122', K122])
    tbl.add_row(['K1111', K1111])
    tbl.add_row(['K1112', K1112])
    tbl.add_row(['K1113', K1113])
    tbl.add_row(['K133', K133])
    print(tbl)
    #print('\n')

def print_constraints(evo, thc, hb, hs, ewpo=None, wmass=None, dm=None):
    tbl = PrettyTable(["Constraint", "Pass/Fail"])
    constraint = {}
    if evo == True:
        constraint['evo'] = 'Pass'
    else:
        constraint['evo'] = 'Fail'
    if thc == True:
        constraint['thc'] = 'Pass'
    else:
        constraint['thc'] = 'Fail'
    if hb == True:
        constraint['hb'] = 'Pass'
    else:
        constraint['hb'] = 'Fail'
    if hs == True:
        constraint['hs'] = 'Pass'
    else:
        constraint['hs'] = 'Fail'
    if ewpo is not None:
        if ewpo == True:
            constraint['ewpo'] = 'Pass'
        else:
            constraint['ewpo'] = 'Fail'
    if wmass is not None:
        if wmass == True:
            constraint['wmass'] = 'Pass'
        else:
            constraint['wmass'] = 'Fail'
    if dm is not None:
        if dm == True:
            constraint['dm'] = 'Pass'
        else:
            constraint['dm'] = 'Fail'
    for key in constraint.keys():
        tbl.add_row([key, constraint[key]])
    print(tbl)

# write the point
def write_valid_point_xsec(runtag, m2, m3, vs, vx, a12, a13, a23, xsec, resfrac):
    outfile = OutputDir + 'trsm_points_' + runtag + '.dat'
    filestream = open(outfile,'a')
    filestream.write(str(m2) + '\t' + str(m3) + '\t' + str(vs) + '\t' + str(vx) + '\t' + str(a12) + '\t' + str(a13) + '\t' + str(a23) + '\t' + str(xsec) + '\t' + str(resfrac) + '\n')
    filestream.close()

def valid_point_info(M2, M3, vs, vx, a12, a13, a23, lX, lPhiX, lSX, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3, evo, thc, hb, hs, ewpo, wmass, dm=None, dm_exclusion_info=None):
    point_info = {
        "M2": M2,
        "M3": M3,
        "vs": vs,
        "vx": vx,
        "a12": a12,
        "a13": a13,
        "a23": a23,
        "lX": lX,
        "lPhiX": lPhiX,
        "lSX": lSX,
        "w1": w1,
        "w2": w2,
        "w3": w3,
        "k1": k1,
        "k2": k2,
        "k3": k3,
        "K111": K111,
        "K112": K112,
        "K113": K113,
        "K123": K123,
        "K122": K122,
        "K1111": K1111,
        "K1112": K1112,
        "K1113": K1113,
        "K133": K133,
        "evo": evo,
        "thc": thc,
        "hb": hb,
        "hs": hs,
        "ewpo": ewpo,
        "wmass": wmass,
        "dm": dm,
    }
    if lPhiX is not None and lSX is not None and vs is not None and a12 is not None:
        point_info["K233"] = lambdas_to_k133_k233(lPhiX, lSX, vs, a12)[1]
    if dm_exclusion_info is not None:
        point_info.update(dm_exclusion_info)
    return point_info


# write the full accepted point record
def write_valid_point(runtag, point_info, MG5xsecs=None):
    outfile = output_path(runtag)
    write_valid_point_file(outfile, point_info, MG5xsecs)

def write_dm_failed_point(runtag, point_info):
    outfile = dm_failed_output_path(runtag)
    write_valid_point_file(outfile, point_info)

def output_path(runtag, suffix=""):
    return OutputDir + 'trsm_points_' + runtag + suffix + '.dat'

def dm_failed_output_path(runtag):
    return output_path(runtag, "_dm_failed")

def reset_output(runtag):
    outfile = output_path(runtag)
    filestream = open(outfile,'w')
    filestream.close()
    if cli_args.write_dm_failed or cli_args.run_ewpt_on_dm_failed:
        filestream = open(dm_failed_output_path(runtag),'w')
        filestream.close()


def ewpt_point_workdir(args, point_index):
    base_dir = args.ewpt_workdir
    if base_dir is None:
        base_dir = Path(OutputDir) / ("ewpt_" + RunTag)
    return Path(base_dir).expanduser() / f"point_{point_index:06d}"


def ewpt_row_from_point_info(point_info):
    return {
        "m1": 125.09,
        "m2": point_info["M2"],
        "m3": point_info["M3"],
        "vs": point_info["vs"],
        "a12": point_info["a12"],
        "lx": point_info["lX"],
        "lphix": point_info["lPhiX"],
        "lsx": point_info["lSX"],
    }


def format_ewpt_candidate(row):
    keys = ("m2", "m3", "vs", "a12", "lx", "lphix", "lsx")
    return "Viable EWPT candidate: " + " ".join(
        f"{key}={row[key]}" for key in keys
    )


def higgstools_analysis_kwargs():
    if not (
        getattr(cli_args, "higgstools_details", False)
        or getattr(cli_args, "save_higgstools_details", True)
    ):
        return {}
    return {
        "print_details": getattr(cli_args, "higgstools_details", False),
        "details_top": getattr(cli_args, "higgstools_top", 5),
        "return_details": getattr(cli_args, "save_higgstools_details", True),
    }


def compact_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def unpack_higgstools_result(result):
    if len(result) >= 3:
        return result[0], result[1], result[2]
    return result[0], result[1], None


def add_higgstools_info(point_info, details):
    if not details:
        return
    higgsbounds = details.get("higgsbounds", {})
    higgssignals = details.get("higgssignals", {})
    point_info["higgstools_hb_selected_limits"] = compact_json(
        higgsbounds.get("selected_limits", {})
    )
    point_info["higgstools_hb_top_obs"] = compact_json(
        higgsbounds.get("top_observed_ratios", [])
    )
    point_info["higgstools_hs_chi2"] = higgssignals.get("chi2")
    point_info["higgstools_hs_delta_chi2"] = higgssignals.get("delta_chi2")
    point_info["higgstools_hs_top_chi2"] = compact_json(
        higgssignals.get("top_chi2_contributors", [])
    )


def effective_m3(args, m2):
    if getattr(args, "resonantDM1", False):
        return 2.0 * float(getattr(args, "m1", 125.09))
    if getattr(args, "resonantDM2", False):
        return 2.0 * float(m2)
    return getattr(args, "m3", None)


def finite_number(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if math.isfinite(number):
        return number
    return None


def one_line_error(error, max_length=1000):
    message = " ".join(str(error).split())
    if len(message) > max_length:
        message = message[:max_length] + "... [truncated]"
    return message


def ewpt_error_details(error):
    sections = [one_line_error(error, max_length=4000)]
    for attribute in ("stdout", "stderr"):
        value = getattr(error, attribute, None)
        if value:
            sections.append(f"\n{attribute}:\n{value}")
    return "\n".join(sections)


def record_ewpt_failure(point_info, workdir, error):
    point_info["ewpt_status"] = "failed"
    point_info["ewpt_error"] = one_line_error(error)
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "ewpt_error.txt").write_text(
        ewpt_error_details(error) + "\n",
        encoding="utf-8",
    )


def select_primary_ewpt_strength(payload):
    strengths = payload.get("transition_strengths") or []
    for temperature_kind in EWPT_STRENGTH_PRIORITY:
        candidates = []
        for strength in strengths:
            if strength.get("temperature_kind") != temperature_kind:
                continue
            ew_true_over_T = finite_number(strength.get("ew_true_over_T"))
            if ew_true_over_T is None:
                continue
            candidates.append((ew_true_over_T, strength))
        if candidates:
            return max(candidates, key=lambda item: item[0])[1]
    return None


def add_ewpt_strength_info(point_info, payload):
    strength = select_primary_ewpt_strength(payload)
    if strength is None:
        return
    point_info["ewpt_ew_true_over_T"] = finite_number(strength.get("ew_true_over_T"))


def add_ewpt_phase_history_info(point_info, payload):
    minimatracer = payload.get("minimatracer") or {}
    global_phase_path = minimatracer.get("global_phase_path") or []
    if isinstance(global_phase_path, str):
        labels = [label.strip() for label in global_phase_path.split("->") if label.strip()]
    else:
        labels = [str(label) for label in global_phase_path]
    if labels:
        point_info["ewpt_global_phase_path"] = " -> ".join(labels)
        point_info["ewpt_has_x_broken"] = any("X_BROKEN" in label for label in labels)
    ew_step_index = minimatracer.get("ew_step_index")
    if ew_step_index is not None:
        point_info["ewpt_ew_step_index"] = ew_step_index


def add_ewpt_info(point_info, payload):
    point_info["ewpt_status"] = "success"
    point_info["ewpt_error"] = ""
    add_ewpt_strength_info(point_info, payload)
    add_ewpt_phase_history_info(point_info, payload)


def run_ewpt_if_requested(
    point_info,
    passed,
    args,
    point_index,
    ewpt_module=None,
    allow_dm_failed=False,
    allow_any_failed=False,
):
    run_for_viable = getattr(args, "run_ewpt", False) and passed
    run_for_dm_failed = (
        getattr(args, "run_ewpt_on_dm_failed", False)
        and allow_dm_failed
        and not passed
    )
    run_for_any_point = getattr(args, "write_all_points", False) and allow_any_failed
    if not (run_for_viable or run_for_dm_failed or run_for_any_point):
        return None

    if ewpt_module is None:
        import test_trsm_ewpt as ewpt_module

    if args.ewpt_require_eq418 and not run_for_any_point:
        ewpt_row = ewpt_row_from_point_info(point_info)
        print(format_ewpt_candidate(ewpt_row))
        eq_4_18 = ewpt_module.check_eq_4_18(ewpt_row)
        if not eq_4_18.satisfied:
            condition_summary = "; ".join(
                f"{name}={value}" for name, value in eq_4_18.conditions.items()
            )
            print(
                "Point is viable but failed Eq. 4.18; skipping EWPT:",
                condition_summary,
            )
            return None

    workdir = ewpt_point_workdir(args, point_index)
    workdir.mkdir(parents=True, exist_ok=True)
    point = ewpt_module.TRSMEWPTPoint(
        index=point_index,
        m2=point_info["M2"],
        m3=point_info["M3"],
        vs=point_info["vs"],
        a12=point_info["a12"],
        lx=point_info["lX"],
        lphix=point_info["lPhiX"],
        lsx=point_info["lSX"],
    )
    config_kwargs = {
        "multistepmode": args.ewpt_multistepmode,
        "thigh": args.ewpt_thigh,
        "plot_phases": args.ewpt_plot_phases,
        "plot_output": args.ewpt_plot_output,
        "plot_format": args.ewpt_plot_format,
        "sym_threshold": args.ewpt_sym_threshold,
        "w1_threshold": args.ewpt_w1_threshold,
        "wx_threshold": args.ewpt_wx_threshold,
        "ws_threshold": args.ewpt_ws_threshold,
    }
    if args.ewpt_executable is not None:
        config_kwargs["executable"] = args.ewpt_executable
    if args.ewpt_minima_executable is not None:
        config_kwargs["minima_executable"] = args.ewpt_minima_executable
    config = ewpt_module.EWPTConfig(**config_kwargs)

    if run_for_any_point and not passed:
        print("Writing all points enabled; running EWPT analysis despite failed constraints in", workdir)
    elif run_for_any_point:
        print("Writing all points enabled; running EWPT analysis in", workdir)
    elif run_for_dm_failed:
        print("Point passes non-DM constraints but fails DM; running EWPT analysis in", workdir)
    else:
        print("Point is viable; running EWPT analysis in", workdir)
    try:
        result = ewpt_module.run_trsm_ewpt(
            point,
            config=config,
            workdir=workdir,
            keep_files=True,
        )
    except Exception as error:
        record_ewpt_failure(point_info, workdir, error)
        print("EWPT analysis failed; continuing scan:", point_info["ewpt_error"])
        print("EWPT error written to", workdir / "ewpt_error.txt")
        return None

    summary = ewpt_module.summarize_result(result)
    payload = ewpt_module.result_to_json(result)
    add_ewpt_info(point_info, payload)
    print(summary)
    (workdir / "ewpt_summary.txt").write_text(summary + "\n", encoding="utf-8")
    (workdir / "ewpt_result.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    print("EWPT summary written to", workdir / "ewpt_summary.txt")
    return result

    
# MAIN FUNCTION:
def evaluate_trsm_point(myseed, m2_val, m3_val, vs_val, vx_val, a12, a13, a23, runmg5=False):
    # get the point information (widths, scalar couplings)
    vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3, h1_BRs, h2_BRs, h3_BRs, xs136_lo_h1, xs136_lo_h2, xs136_lo_h3 = generate_lams(myseed, m2_val, m3_val, vs_val, vx_val, a12, a13, a23, PRINTINFO)
    Lambdas =[K111,K112,K113,K123,K122,K1111,K1112,K1113,K133]
    if debug is True:
        print_info(vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3)

    # check EWPO
    sinth = np.sin(a12)
    #EWPO_cur_old = check_EWPO(125.09, M2, sinth, Mz, Mw, Delta_S_central, Delta_T_central, errS, errT, covST) #current # U=0
    EWPO_cur = check_EWPO_wU(125.09, M2, sinth, Mz, Mw, Delta_S_central_wU, Delta_T_central_wU, Delta_U_central_wU, errS_wU, errT_wU, errU_wU, covST_wU, covSU_wU, covTU_wU) #current # U!=0

    #if EWPO_cur != EWPO_cur_wU:
    #    print("EWPOWARNING", EWPO_cur, EWPO_cur_wU)
    

    # check W mass:
    wmass = check_wmass_tania(M2, sinth)
    
        
    # check HiggsTools:
    hb, hs, higgstools_details = unpack_higgstools_result(analyze_parampoint(
        pred,
        H1,
        H2,
        H3,
        125.09,
        M2,
        M3,
        k1,
        k2,
        k3,
        h1_BRs,
        h2_BRs,
        h3_BRs,
        **higgstools_analysis_kwargs(),
    ))
    if debug is False:
        if hb is False or hs is False:
            return 0
    #thc = theory_constraints(vs, vx, M2, M3, a12, a13, a23)
    #if debug is False:
    #    if thc is False:
    #        return 0
    # test the cosmological constraints
    #evo = test_evo(vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133)
    #if debug is False:
    #    if evo is False:
    #        return 0

    # temporarily removing theory constraints: to be REINSTATED!
    evo = True
    thc = True
    if debug is True:
        print_constraints(evo, thc, hb, hs, EWPO_cur, wmass)
    # get the hh cross section
    # if all constraints are ok, check the xsec for hhh:
    if evo is True and thc is True and hb is True and hs is True and EWPO_cur is True and wmass is True:
        if debug is False:
            print_info(vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3)
            print_constraints(evo, thc, hb, hs, EWPO_cur, wmass)
        MG5xsecs = {}
        if runmg5 is True:
            print('All constraints passed, running selected MG5 processes, please wait!')
            MG5xsecs = run_mg5_processes(MG5ProcessesToRun, 'SCAN' + str(Energy), Lambdas, k1, k2, k3, M2, w2, M3, w3, Energy)
            print('MG5 cross sections [pb] =', MG5xsecs)
        point_info = valid_point_info(M2, M3, vs, vx, a12, a13, a23, None, None, None, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3, evo, thc, hb, hs, EWPO_cur, wmass)
        add_higgstools_info(point_info, higgstools_details)
        write_valid_point(RunTag, point_info, MG5xsecs)
        return 1
    return 0

# MAIN FUNCTION for vx=0:
def evaluate_trsm_point_vxzero(myseed, m2_val, m3_val, vs_val, a12, lX, lPhiX, lSX, runmg5=False, report=False, write_all=False, point_index=1):
    write_all_points = write_all or getattr(cli_args, "write_all_points", False)
    force_ewpt_for_all_points = getattr(cli_args, "write_all_points", False)
    print_info_enabled = getattr(cli_args, "print_info", False)
    short_circuit_failures = (
        debug is False
        and report is False
        and write_all_points is False
        and print_info_enabled is False
    )
    # fix a23, a13, vx to zero:
    vx_val = 0
    a13 = 0
    a23 = 0
    # get the point information (widths, scalar couplings)
    vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3, h1_BRs, h2_BRs, h3_BRs, xs136_lo_h1, xs136_lo_h2, xs136_lo_h3 = generate_lams(myseed, m2_val, m3_val, vs_val, vx_val, a12, a13, a23, PRINTINFO, lX=lX, lPhiX=lPhiX, lSX=lSX)
    Lambdas =[K111,K112,K113,K123,K122,K1111,K1112,K1113,K133]
    if debug is True or report is True or print_info_enabled is True:
        print_info_vxzero(vs, vx, M2, M3, a12, a13, a23, lX, lPhiX, lSX, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3)

    # check EWPO
    sinth = np.sin(a12)
    #EWPO_cur_old = check_EWPO(125.09, M2, sinth, Mz, Mw, Delta_S_central, Delta_T_central, errS, errT, covST) #current # U=0
    EWPO_cur = check_EWPO_wU(125.09, M2, sinth, Mz, Mw, Delta_S_central_wU, Delta_T_central_wU, Delta_U_central_wU, errS_wU, errT_wU, errU_wU, covST_wU, covSU_wU, covTU_wU) #current # U!=0

    #if EWPO_cur != EWPO_cur_wU:
    #    print("EWPOWARNING", EWPO_cur, EWPO_cur_wU)
    

    # check W mass:
    wmass = check_wmass_tania(M2, sinth)
    
    # check HiggsTools:
    hb, hs, higgstools_details = unpack_higgstools_result(analyze_parampoint(
        pred,
        H1,
        H2,
        H3,
        125.09,
        M2,
        M3,
        k1,
        k2,
        k3,
        h1_BRs,
        h2_BRs,
        h3_BRs,
        **higgstools_analysis_kwargs(),
    ))
    if short_circuit_failures:
        if hb is False or hs is False:
            return 0
    thc = theory_constraints_vxzero(vs, M2, M3, a12, lX, lPhiX, lSX)
    if short_circuit_failures:
        if thc is False:
            return 0
    # test the cosmological constraints
    evo = test_evo_vxzero(vs, M2, M3, a12, lX, lPhiX, lSX, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133)
    if short_circuit_failures:
        if evo is False:
            return 0
    dm = test_dm(lX, lPhiX, lSX, M3, vs, a12, M2)
    if debug is True or report is True or print_info_enabled is True:
        print_constraints(evo, thc, hb, hs, EWPO_cur, wmass, dm[0])
        print_dm_info(dm[1])
    pre_dm_passed = evo is True and thc is True and hb is True and hs is True and EWPO_cur is True and wmass is True
    point_info = valid_point_info(M2, M3, vs, vx, a12, a13, a23, lX, lPhiX, lSX, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3, evo, thc, hb, hs, EWPO_cur, wmass, dm[0], dm[2])
    add_higgstools_info(point_info, higgstools_details)
    dm_failed_but_otherwise_allowed = pre_dm_passed and dm[0] is False
    if dm_failed_but_otherwise_allowed:
        ewpt_error = None
        if cli_args.run_ewpt_on_dm_failed and force_ewpt_for_all_points is False:
            try:
                run_ewpt_if_requested(
                    point_info,
                    False,
                    cli_args,
                    point_index,
                    allow_dm_failed=True,
                )
            except Exception as error:
                ewpt_error = error
        if cli_args.write_dm_failed or cli_args.run_ewpt_on_dm_failed:
            write_dm_failed_point(RunTag, point_info)
            print('Point passes non-DM constraints but fails DM; written to', dm_failed_output_path(RunTag))
        if ewpt_error is not None:
            raise ewpt_error
    if short_circuit_failures:
        if dm[0] is False:
            return 0
    # get the hh cross section
    # if all constraints are ok, check the xsec for hhh:
    passed = pre_dm_passed and dm[0] is True
    if passed or write_all_points is True:
        if debug is False and report is False and print_info_enabled is False:
            print_info_vxzero(vs, vx, M2, M3, a12, a13, a23, lX, lPhiX, lSX, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3)
            print_constraints(evo, thc, hb, hs, EWPO_cur, wmass, dm[0])
            print_dm_info(dm[1])
        MG5xsecs = {}
        if passed and runmg5 is True:
            print('All constraints passed, running selected MG5 processes, please wait!')
            MG5xsecs = run_mg5_processes(MG5ProcessesToRun, 'SCAN' + str(Energy), Lambdas, k1, k2, k3, M2, w2, M3, w3, Energy)
            print('MG5 cross sections [pb] =', MG5xsecs)
        ewpt_error = None
        if passed or force_ewpt_for_all_points:
            try:
                ewpt_kwargs = {}
                if force_ewpt_for_all_points:
                    ewpt_kwargs["allow_any_failed"] = True
                run_ewpt_if_requested(
                    point_info,
                    passed,
                    cli_args,
                    point_index,
                    **ewpt_kwargs,
                )
            except Exception as error:
                ewpt_error = error
        write_valid_point(RunTag, point_info, MG5xsecs)
        if report is True:
            print('Point record written to', output_path(RunTag))
        if ewpt_error is not None:
            raise ewpt_error
    return 1 if passed else 0

# round to sgf significant figures
def round_signif(m2, m3, vs, vx, a12, a13, a23, lX, lPhiX, lSX, sgf):
    return round_sig(m2, sgf), round_sig(m3, sgf), round_sig(vs, sgf), 0, round_sig(a12, sgf), 0, 0, round_sig(lX, sgf), round_sig(lPhiX, sgf), round_sig(lSX, sgf) 


# random number either -1 or 1:
def randsign():
    return 1 if random.random() < 0.5 else -1


############################################################
# define ranges here
############################################################

# ranges of masses m2 and m3 to scan over
m2_min=4
m2_max=20

m3_min=500
m3_max=700

# ranges of vevs
vs_min=400
vs_max=1000

# no range for vx (DM candidate)!
vx_min=0
vx_max=0

# ranges of k1, k2, k3 (can be positive or negative)
k1_min = 0.95
k1_max = 1.0

# for the grid scan:
num_m2 = 2
num_m3 = 2

# ranges of couplings if vx=0
lX_min = 0.0
lX_max= 0.5

lPhiX_min = 0.0
lPhiX_max = 0.002

lSX_min = 0.0
lSX_max = 0.002

# ranges of physical dimensionful couplings if --scan-k133-k233 is used [GeV]
K133_min = -1.0
K133_max = 1.0

K233_min = -1.0
K233_max = 1.0


############################################################
# Scan begins here
############################################################

def has_explicit_point(args):
    if args.resonantDM1 or args.resonantDM2:
        point_args = [name for name in TRSM_POINT_ARGS if name != "m3"]
    else:
        point_args = TRSM_POINT_ARGS
    return all(getattr(args, name) is not None for name in point_args)


def run_single_vxzero_point(args):
    m3 = effective_m3(args, args.m2)
    print('Evaluating one explicit vx=0 TRSM point')
    print(
        'm2 =', args.m2,
        'm3 =', m3,
        'vs =', args.vs,
        'a12 =', args.a12,
        'lx =', args.lx,
        'lphix =', args.lphix,
        'lsx =', args.lsx,
    )
    evalpoint = evaluate_trsm_point_vxzero(
        ini_seed,
        args.m2,
        m3,
        args.vs,
        args.a12,
        args.lx,
        args.lphix,
        args.lsx,
        runmg5=RunMG5,
        report=True,
        write_all=True,
        point_index=1,
    )
    print('Evaluated 1 explicit point, out of which', evalpoint, 'are viable')
    return evalpoint


def run_random_vxzero_scan():
    # SCAN:
    passcounter = 0 # count number of passing points
    print('Generating points randomly within ranges')
    random.seed(ini_seed)

    # fixed parameters:
    vx = 0
    k3=0
    a13 = 0
    a23 = 0

    for i in tqdm(range(0,nrandom)):
        # scan over free parameters
        k1=random.uniform(k1_min,k1_max)
        m2=random.uniform(m2_min, m2_max)
        if cli_args.resonantDM1 or cli_args.resonantDM2:
            m3=effective_m3(cli_args, m2)
        else:
            m3=random.uniform(m2+mhiggs, m3_max)
        vs=random.uniform(vs_min, vs_max)
        # free parameters for vx=0:
        lX=random.uniform(lX_min, lX_max)
        if getattr(cli_args, "scan_k133_k233", False):
            K133=random.uniform(K133_min, K133_max)
            K233=random.uniform(K233_min, K233_max)
            lPhiX=0
            lSX=0
        else:
            lPhiX=random.uniform(lPhiX_min, lPhiX_max)
            lSX=random.uniform(lSX_min, lSX_max)

        # dependent parameters
        a12 = randsign() * np.arccos(k1)
        k2=randsign() * np.sqrt(1-k1**2)

        # optional (vx!=0), convert angles to k1, k2, k3:
        #a12, a13, a23 = ks_to_angles(k1,k2,k3)

        # round to 4 significant figures:
        m2, m3, vs, vx, a12, a13, a23, lX, lPhiX, lSX = round_signif(m2, m3, vs, vx, a12, a13, a23, lX, lPhiX, lSX, 4)
        if getattr(cli_args, "scan_k133_k233", False):
            K133=round_sig(K133, 4)
            K233=round_sig(K233, 4)
            lPhiX, lSX = k133_k233_to_lambdas(K133, K233, vs, a12)
        if cli_args.resonantDM1 or cli_args.resonantDM2:
            m3=round_sig(effective_m3(cli_args, m2), 4)

        # evaluate: evalpoint is 1 if point passes, 0 if not
        if vx != 0:
                evalpoint = evaluate_trsm_point(ini_seed, m2, m3, vs, vx, a12, a13, a23,runmg5=RunMG5)
        else:
                # if point passes, file will be written with order: m2, m3, vs, a12, lX, lPhiX, lSX + any additional info to be determined
                evalpoint = evaluate_trsm_point_vxzero(ini_seed, m2, m3, vs, a12, lX, lPhiX, lSX,runmg5=RunMG5, point_index=i + 1)
        # count the passing points:
        passcounter = passcounter + evalpoint

    print('Generated', nrandom,'points, out of which', passcounter, 'are viable')
    return passcounter


def main():
    global RunTag

    print('\nScanning TRSM parameter space')
    RunTag = RunTag + '_vxzero'
    if has_explicit_point(cli_args):
        RunTag = RunTag + '_manual'

    # reset the output file?
    if ResetOutput is True:
        reset_output(RunTag)

    if has_explicit_point(cli_args):
        return run_single_vxzero_point(cli_args)
    return run_random_vxzero_scan()


if __name__ == "__main__":
    main()
