import math
from math import log10, floor
from pathlib import Path
import Higgs.predictions as HP
import Higgs.bounds as HB
import Higgs.signals as HS
import numpy as np
import pandas as pd
from prettytable import PrettyTable
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
from tqdm import tqdm
from matplotlib import ticker, cm
from scipy import interpolate
from scipy.optimize import fsolve, root
from functools import partial

def round_sig(x, sig=2):
    if x == 0.:
        return 0.
    if math.isnan(x) is True:
        print('Warning, NaN!')
        return 0.
    return round(x, sig-int(floor(log10(abs(x))))-1)


##################################
# SET UP HIGGS TOOLS:
#####################################

mhSM = 125.09
PROJECT_ROOT = Path(__file__).resolve().parents[1]

pred = HP.Predictions() # create the model predictions

bounds = HB.Bounds(str(PROJECT_ROOT / 'hbdataset')) # load HB dataset
signals = HS.Signals(str(PROJECT_ROOT / 'hsdataset')) # load HS dataset

# add a SM-like Higgs boson with SM-like couplings

H1 = pred.addParticle(HP.BsmParticle("H1", "neutral", "even"))

H1.setMass(mhSM)
HP.effectiveCouplingInput(H1, HP.scaledSMlikeEffCouplings(1.0))
# get the SM chi-squared for HiggsSignals
ress_SM = signals(pred)

print("HiggsSignals chi-sq. for SM=", ress_SM)

# add second BSM Higgs boson 

H2 = pred.addParticle(HP.BsmParticle("H2", "neutral", "even"))

# add third BSM Higgs boson

H3 = pred.addParticle(HP.BsmParticle("H3", "neutral", "even"))


##########################################
# FUNCTIONS 
##########################################

# ensure that the BR arrays sum up to unity to avoid issues with HiggsTools
def ensure_sum_unit(brarray):
    sumarray = brarray.sum() - brarray[-1]
    brarray_new = np.array([None] * len(brarray))
    for i in range(len(brarray)-1):
        brarray_new[i] = round_sig(brarray[i]/sumarray,100)
    brarray_new[-1] = round_sig(brarray[-1],3)
    #print(brarray_new.sum() - brarray_new[-1])
    return brarray_new


def _method_value(obj, name, default=""):
    try:
        value = getattr(obj, name)
    except AttributeError:
        return default
    if callable(value):
        try:
            value = value()
        except TypeError:
            return default
    return value


def _display_value(value):
    if value is None:
        return ""
    if hasattr(value, "name"):
        return value.name
    return str(value)


def _format_float(value):
    try:
        return f"{float(value):.6g}"
    except (TypeError, ValueError):
        return ""


def _float_or_none(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if math.isfinite(number):
        return number
    return None


def _int_or_string(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return _display_value(value)


def applied_limit_summary(applied_limit):
    limit = _method_value(applied_limit, "limit", None)
    particles = _method_value(applied_limit, "contributingParticles", [])
    if not isinstance(particles, (list, tuple)):
        particles = [particles] if particles else []
    return {
        "obsRatio": _float_or_none(_method_value(applied_limit, "obsRatio")),
        "expRatio": _float_or_none(_method_value(applied_limit, "expRatio")),
        "particles": [str(particle) for particle in particles],
        "limit_id": _int_or_string(_method_value(limit, "id"))
        if limit is not None
        else "",
        "reference": _method_value(limit, "reference") if limit is not None else "",
        "cite_key": _method_value(limit, "citeKey") if limit is not None else "",
        "process": _method_value(limit, "processDesc") if limit is not None else "",
    }


def _applied_limit_row(applied_limit, label=""):
    summary = applied_limit_summary(applied_limit)
    return [
        label,
        _format_float(summary["obsRatio"]),
        _format_float(summary["expRatio"]),
        ",".join(summary["particles"]),
        summary["limit_id"],
        summary["reference"],
        summary["cite_key"],
        summary["process"],
    ]


def higgsbounds_summary(result, top_n=5):
    selected_limits = getattr(result, "selectedLimits", {}) or {}
    applied_limits = list(getattr(result, "appliedLimits", []) or [])
    applied_limits.sort(
        key=lambda item: _method_value(item, "obsRatio", float("-inf")),
        reverse=True,
    )
    return {
        "allowed": bool(result.allowed),
        "selected_limits": {
            str(particle): applied_limit_summary(applied_limit)
            for particle, applied_limit in selected_limits.items()
        },
        "top_observed_ratios": [
            applied_limit_summary(applied_limit)
            for applied_limit in applied_limits[:top_n]
        ],
    }


def print_higgsbounds_details(result, top_n=5):
    print("HiggsBounds allowed =", result.allowed)

    selected_limits = getattr(result, "selectedLimits", {}) or {}
    if selected_limits:
        table = PrettyTable(
            [
                "particle",
                "obsRatio",
                "expRatio",
                "for",
                "limit id",
                "reference",
                "cite key",
                "process",
            ]
        )
        for particle, applied_limit in selected_limits.items():
            table.add_row(_applied_limit_row(applied_limit, label=particle))
        print("HiggsBounds selected limits")
        print(table)

    applied_limits = list(getattr(result, "appliedLimits", []) or [])
    applied_limits.sort(
        key=lambda item: _method_value(item, "obsRatio", float("-inf")),
        reverse=True,
    )
    if applied_limits:
        table = PrettyTable(
            [
                "rank",
                "obsRatio",
                "expRatio",
                "for",
                "limit id",
                "reference",
                "cite key",
                "process",
            ]
        )
        for rank, applied_limit in enumerate(applied_limits[:top_n], start=1):
            table.add_row(_applied_limit_row(applied_limit, label=rank))
        print(f"HiggsBounds top {min(top_n, len(applied_limits))} observed ratios")
        print(table)


def higgs_signals_chisq_rows(signals_obj, pred):
    rows = []
    for measurement in signals_obj.measurements():
        try:
            contributions = measurement.chisqContributions(pred)
        except Exception:
            continue
        submeasurements = _method_value(measurement, "subMeasurements", {}) or {}
        for sub_name, chisq in contributions.items():
            submeasurement = submeasurements.get(sub_name)
            rows.append(
                {
                    "chisq": float(chisq),
                    "measurement_id": _method_value(measurement, "id"),
                    "reference": _method_value(measurement, "reference"),
                    "cite_key": _method_value(measurement, "citeKey"),
                    "experiment": _display_value(_method_value(measurement, "experiment")),
                    "collider": _display_value(_method_value(measurement, "collider")),
                    "luminosity": _method_value(measurement, "luminosity"),
                    "submeasurement": sub_name,
                    "process": _method_value(submeasurement, "processDesc")
                    if submeasurement is not None
                    else "",
                }
            )
    rows.sort(key=lambda row: row["chisq"], reverse=True)
    return rows


def higgssignals_summary(pred, chi2, top_n=5):
    rows = higgs_signals_chisq_rows(signals, pred)
    return {
        "chi2": _float_or_none(chi2),
        "sm_chi2": _float_or_none(ress_SM),
        "delta_chi2": _float_or_none(chi2 - ress_SM),
        "top_chi2_contributors": rows[:top_n],
    }


def print_higgssignals_details(pred, chi2, top_n=5):
    delta_chi2 = chi2 - ress_SM
    print("HiggsSignals chi2 =", chi2)
    print("HiggsSignals SM chi2 =", ress_SM)
    print("HiggsSignals delta chi2 =", delta_chi2)

    rows = higgs_signals_chisq_rows(signals, pred)
    if not rows:
        return
    table = PrettyTable(
        [
            "rank",
            "chi2",
            "measurement id",
            "reference",
            "cite key",
            "experiment",
            "collider",
            "lumi",
            "submeasurement",
            "process",
        ]
    )
    for rank, row in enumerate(rows[:top_n], start=1):
        table.add_row(
            [
                rank,
                _format_float(row["chisq"]),
                row["measurement_id"],
                row["reference"],
                row["cite_key"],
                row["experiment"],
                row["collider"],
                _format_float(row["luminosity"]),
                row["submeasurement"],
                row["process"],
            ]
        )
    print(f"HiggsSignals top {min(top_n, len(rows))} chi2 contributors")
    print(table)


def print_higgstools_details(pred, bounds_result, signals_chi2, top_n=5):
    print_higgsbounds_details(bounds_result, top_n=top_n)
    print_higgssignals_details(pred, signals_chi2, top_n=top_n)


def higgstools_summary(pred, bounds_result, signals_chi2, top_n=5):
    return {
        "higgsbounds": higgsbounds_summary(bounds_result, top_n=top_n),
        "higgssignals": higgssignals_summary(pred, signals_chi2, top_n=top_n),
    }


def analyze_parampoint(
    pred,
    H1,
    H2,
    H3,
    M1,
    M2,
    M3,
    k1,
    k2,
    k3,
    h1_BRs,
    h2_BRs,
    h3_BRs,
    print_details=False,
    details_top=5,
    return_details=False,
):

    h1_BRs = ensure_sum_unit(h1_BRs)
    h2_BRs = ensure_sum_unit(h2_BRs)

    # set the present masses for the scalars
    H1.setMass(M1)
    H2.setMass(M2)
    # get the branching ratios
    HP.effectiveCouplingInput(H1, HP.scaledSMlikeEffCouplings(k1))
    HP.effectiveCouplingInput(H2, HP.scaledSMlikeEffCouplings(k2))
    # RESET BRs BEFORE SETTING THEM TO AVOID ISSUES WITH BR>1
    # H -> bbbar | H -> tautau | H -> mumu | H -> cc | H -> ss | H -> tt | H -> gg | H -> gammagamma | H -> Zgamma | H -> WW | H -> ZZ
    # H1
    H1.setBr('bb', 0.)
    H1.setBr('tautau', 0.)
    H1.setBr('mumu', 0.)
    H1.setBr('cc', 0.)
    H1.setBr('ss', 0.)
    H1.setBr('tt', 0.)
    H1.setBr('gg', 0.)
    H1.setBr('gamgam', 0.)
    H1.setBr('Zgam', 0.)
    H1.setBr('WW', 0.)
    H1.setBr('ZZ', 0.)
    # SET THE BRS
    H1.setBr('bb', h1_BRs[0])
    H1.setBr('tautau', h1_BRs[1])
    H1.setBr('mumu', h1_BRs[2])
    H1.setBr('cc', h1_BRs[3])
    H1.setBr('ss', h1_BRs[4])
    H1.setBr('tt', h1_BRs[5])
    H1.setBr('gg', h1_BRs[6])
    H1.setBr('gamgam', h1_BRs[7])
    H1.setBr('Zgam',h1_BRs[8])
    H1.setBr('WW', h1_BRs[9])
    H1.setBr('ZZ', h1_BRs[10])
    H1.setTotalWidth(h1_BRs[-1])
    # H2
    # RESET BRs BEFORE SETTING THEM TO AVOID ISSUES WITH BR>1
    H2.setBr('bb', 0.)
    H2.setBr('tautau', 0.)
    H2.setBr('mumu', 0.)
    H2.setBr('cc', 0.)
    H2.setBr('ss', 0.)
    H2.setBr('tt', 0.)
    H2.setBr('gg', 0.)
    H2.setBr('gamgam', 0.)
    H2.setBr('Zgam', 0.)
    H2.setBr('WW', 0.)
    H2.setBr('ZZ', 0.)
    # SET THE BRS
    H2.setBr('bb', h2_BRs[0])
    H2.setBr('tautau', h2_BRs[1])
    H2.setBr('mumu', h2_BRs[2])
    H2.setBr('cc', h2_BRs[3])
    H2.setBr('ss', h2_BRs[4])
    H2.setBr('tt', h2_BRs[5])
    H2.setBr('gg', h2_BRs[6])
    H2.setBr('gamgam', h2_BRs[7])
    H2.setBr('Zgam',h2_BRs[8])
    H2.setBr('WW',h2_BRs[9])
    H2.setBr('ZZ',h2_BRs[10])
    H2.setBr('H1','H1',h2_BRs[11])
    H2.setTotalWidth(h2_BRs[-1])

    # get the HiggsBounds result
    resb = bounds(pred)

    # get the HiggsSignals result
    ress = signals(pred)
    #print(signals(pred).appliedLimits)
    #print(ress)
    if print_details:
        print_higgstools_details(pred, resb, ress, top_n=details_top)

    HS_allowed = False
    if ress - ress_SM < 4.00:
        HS_allowed = True
    else:
        HS_allowed = False

    if return_details:
        details = higgstools_summary(pred, resb, ress, top_n=details_top)
        return resb.allowed, HS_allowed, details
    return resb.allowed, HS_allowed
