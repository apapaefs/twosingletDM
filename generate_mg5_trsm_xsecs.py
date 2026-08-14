import fcntl
import gzip
import math
import os
import subprocess
from contextlib import contextmanager
from math import floor, log10
from pathlib import Path


# MG5/aMC subdirectory.  Override this for installations outside the project.
MGLocation = os.environ.get(
    "TRSM_MG5_LOCATION",
    str(Path(__file__).resolve().parents[1] / "MG5_aMC_v3_5_15") + os.sep,
)


# Generated process directories, relative to MGLocation.
ProcLocation = {
    "hh": "gg_hh_twoscalar/",
    "hhh": "gg_hhh_twoscalar/",
    "gg_heta0": "gg_heta0/",
    "pp_eta0Z": "pp_eta0Z/",
}


TREE_LEVEL_SURVEY_ITERATIONS = 3
LOOP_INDUCED_SURVEY_ITERATIONS = 1


def round_sig(x, sig=2):
    if x == 0.0:
        return 0.0
    if math.isnan(x):
        print("Warning, NaN!", x)
        return 0.0
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


def _process_directory(process, mgloc):
    if process not in ProcLocation:
        available = ", ".join(sorted(ProcLocation))
        raise ValueError(
            f"MG5 process {process!r} is not defined; available processes: {available}"
        )
    process_dir = Path(mgloc).expanduser() / ProcLocation[process]
    madevent = process_dir / "bin" / "madevent"
    if not process_dir.is_dir() or not madevent.is_file():
        raise FileNotFoundError(
            f"Generated MG5 process directory for {process!r} was not found at "
            f"{process_dir}. Generate/copy that process there or set "
            "TRSM_MG5_LOCATION to the MG5 installation containing it."
        )
    return process_dir, madevent


def _run_name(runnum, m2, w2, m3, w3, lambdas):
    return (
        "run"
        + str(runnum)
        + "_m2_"
        + str(m2)
        + "_m3_"
        + str(m3)
        + "_w2_"
        + str(w2)
        + "_w3_"
        + str(w3)
        + "_"
        + "_".join(lambdas)
    )


def _lhe_path(process_dir, run_name):
    return process_dir / "Events" / run_name / "unweighted_events.lhe.gz"


def _survey_iterations(process_dir):
    """Use MadEvent's minimum survey length for tree-level processes."""
    characteristics = process_dir / "SubProcesses" / "proc_characteristics"
    try:
        lines = characteristics.read_text(encoding="ascii").splitlines()
    except OSError:
        return TREE_LEVEL_SURVEY_ITERATIONS

    for line in lines:
        key, separator, value = line.partition("=")
        if separator and key.strip() == "loop_induced":
            if value.strip().lower() == "true":
                return LOOP_INDUCED_SURVEY_ITERATIONS
            return TREE_LEVEL_SURVEY_ITERATIONS
    return TREE_LEVEL_SURVEY_ITERATIONS


@contextmanager
def _madevent_lock(process_dir):
    """Serialize runs that share MadEvent's mutable process directory."""
    lock_path = process_dir / ".trsm_madevent.lock"
    with lock_path.open("a+", encoding="ascii") as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def drive_mg(
    process,
    runnum,
    mgloc,
    k1choice,
    k2choice,
    k3choice,
    LambdasArray,
    m2,
    w2,
    m3,
    w3,
    nevents,
    nruns,
    output=False,
    ecm=13,
    *,
    w1=None,
    k233=None,
):
    """Run a generated MadEvent process for the supplied TRSM parameter card."""
    process_dir, madevent = _process_directory(process, mgloc)
    command_file = process_dir / f"mg5_{process}_lambdavar_run{runnum}.dcmd"
    ebeam = ecm * 1000.0 / 2.0
    survey_iterations = _survey_iterations(process_dir)
    counter = 0
    with _madevent_lock(process_dir):
        for lambdas in LambdasArray:
            if counter >= nruns:
                break
            run_name = _run_name(runnum, m2, w2, m3, w3, lambdas)
            lhefile = _lhe_path(process_dir, run_name)
            if lhefile.exists():
                counter += 1
                continue

            commands = [
                f"generate_events {run_name} --accuracy=0.25 --points=300 "
                f"--iterations={survey_iterations}",
                f"set ebeam1 {ebeam}",
                f"set ebeam2 {ebeam}",
                f"set Meta {m2}",
                f"set Weta {w2}",
                f"set Miota {m3}",
                f"set Wiota {w3}",
                f"set k1 {k1choice}",
                f"set k2 {k2choice}",
                f"set k3 {k3choice}",
                f"set kap111 {lambdas[0]}",
                f"set kap112 {lambdas[1]}",
                f"set kap113 {lambdas[2]}",
                f"set kap123 {lambdas[3]}",
                f"set kap122 {lambdas[4]}",
                f"set kap1111 {lambdas[5]}",
                f"set kap1112 {lambdas[6]}",
                f"set kap1113 {lambdas[7]}",
                f"set kap133 {lambdas[8]}",
            ]
            if w1 is not None:
                commands.append(f"set WH {w1}")
            effective_k233 = k233 if k233 is not None else (
                lambdas[9] if len(lambdas) > 9 else None
            )
            if effective_k233 is not None:
                commands.append(f"set kap233 {effective_k233}")
            commands.extend((f"set nevents {nevents}", "0"))
            command_file.write_text("\n".join(commands), encoding="ascii")

            if output:
                print(command_file.read_text(encoding="ascii"))
            completed = subprocess.run(
                [str(madevent), str(command_file)],
                cwd=process_dir,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            if output and completed.stdout:
                print(completed.stdout)
            if completed.returncode != 0 or not lhefile.exists():
                tail = "\n".join((completed.stdout or "").splitlines()[-30:])
                raise RuntimeError(
                    f"MG5 process {process!r} failed for run {run_name} "
                    f"(exit code {completed.returncode}); expected {lhefile}.\n{tail}"
                )
            counter += 1
    print("Done generating cross section")
    return counter


def _integrated_weight(lhefile):
    with gzip.open(lhefile, "rt", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if "Integrated weight" not in line:
                continue
            value_text = line.rsplit(":", 1)[-1].strip().split()[0]
            return float(value_text)
    raise ValueError(f"No Integrated weight entry was found in {lhefile}")


def read_files(runnum, LambdasArray, m2, w2, m3, w3, process, nruns):
    process_dir, _madevent = _process_directory(process, MGLocation)
    X = []
    Z = []
    XSEC = {}
    counter = 0
    for lambdas in LambdasArray:
        if counter >= nruns:
            break
        run_name = _run_name(runnum, m2, w2, m3, w3, lambdas)
        lhefile = _lhe_path(process_dir, run_name)
        if not lhefile.exists():
            raise FileNotFoundError(f"MG5 LHE file does not exist: {lhefile}")
        xsec = _integrated_weight(lhefile)
        lambdas_tuple = tuple(float(value) for value in lambdas)
        X.append(lambdas_tuple)
        Z.append(xsec)
        XSEC[lambdas_tuple] = xsec
        counter += 1
    return X, Z, XSEC


def get_mg5_xsec(
    process,
    runnum,
    LambdasArray,
    k1,
    k2,
    k3,
    m2,
    w2,
    m3,
    w3,
    ecm=13,
    *,
    w1=None,
    k233=None,
):
    rounded_lambdas = [str(round_sig(value, 4)) for value in LambdasArray]
    rounded_m2 = round_sig(m2, 4)
    rounded_w2 = round_sig(w2, 4)
    rounded_m3 = round_sig(m3, 4)
    rounded_w3 = round_sig(w3, 4)
    drive_mg(
        process,
        runnum,
        MGLocation,
        round_sig(k1, 4),
        round_sig(k2, 4),
        round_sig(k3, 4),
        [rounded_lambdas],
        rounded_m2,
        rounded_w2,
        rounded_m3,
        rounded_w3,
        1,
        1,
        output=True,
        ecm=ecm,
        w1=round_sig(w1, 4) if w1 is not None else None,
        k233=round_sig(k233, 4) if k233 is not None else None,
    )
    _x, cross_sections, _mapping = read_files(
        runnum,
        [rounded_lambdas],
        rounded_m2,
        rounded_w2,
        rounded_m3,
        rounded_w3,
        process,
        1,
    )
    if not cross_sections:
        raise RuntimeError(f"MG5 returned no cross section for process {process!r}")
    return cross_sections[0]
