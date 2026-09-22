"""Planck 2018 energy-injection constraint in the micrOMEGAs s-wave approximation."""

import json
import math
import subprocess
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path


CMB_METHOD = "micromegas_planck2018_swave_v1"
CMB_OUTPUT_PREFIX = "TRSM_PlanckCMB_v1 "
CMB_PANN_LIMIT = 3.2e-28  # cm^3 s^-1 GeV^-1, 95% CL
CMB_OMEGA_REFERENCE = 0.12
CMB_COLUMNS = (
    "dm_cmb_enabled", "dm_cmb_available", "dm_cmb_status", "dm_cmb_reason",
    "dm_cmb_ratio_raw", "dm_cmb_abundance_fraction", "dm_cmb_ratio",
    "dm_cmb_excluded",
)


def add_cmb_arguments(parser):
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--planck-cmb", dest="planck_cmb", action="store_true",
        help="Apply the relic-rescaled Planck CMB annihilation bound (default for new version 7 runs).",
    )
    group.add_argument(
        "--no-planck-cmb", dest="planck_cmb", action="store_false",
        help="Disable the Planck CMB bound (default for version 6 and legacy resumes).",
    )
    parser.set_defaults(planck_cmb=None)


def cmb_configuration(args):
    enabled = bool(getattr(args, "planck_cmb", False))
    if not enabled:
        return {"enabled": False}
    result = {
        "enabled": True, "method": CMB_METHOD,
        "pann_limit_cm3_s_GeV": CMB_PANN_LIMIT, "confidence_level": 0.95,
        "omega_reference": CMB_OMEGA_REFERENCE,
        "abundance_rescaling": "min(1, Omega_h2 / 0.12)^2",
        "approximation": "built-in low-velocity s-wave; calcSpectrum default V0=sqrt(3)*vRot/c",
        "source": "https://arxiv.org/html/2606.06645v1",
    }
    driver = getattr(args, "_cmb_driver", None)
    if driver is not None:
        result["driver"] = driver
    return result


@lru_cache(maxsize=16)
def _driver_capability(path, mtime_ns, size):
    try:
        completed = subprocess.run(
            [path, "--capabilities"], capture_output=True, text=True,
            check=True, timeout=10,
        )
        payload = json.loads(completed.stdout)
        capability = payload["planck_cmb"]
        if payload.get("schema") != "trsm_driver_capabilities_v1":
            raise ValueError("unsupported capability schema")
        if capability.get("method") != CMB_METHOD:
            raise ValueError("unsupported CMB method")
        if capability.get("pann_limit_cm3_s_GeV") != CMB_PANN_LIMIT:
            raise ValueError("unexpected CMB normalization")
        if capability.get("spectrum_key") != 7:
            raise ValueError("unexpected spectrum calculation")
        if capability.get("vz_decay") != 0 or capability.get("vw_decay") != 0:
            raise ValueError("unexpected off-shell spectrum settings")
        for key in ("spectra_flag", "vrot_km_s"):
            value = capability[key]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"invalid {key}")
        if capability["vrot_km_s"] <= 0:
            raise ValueError("invalid annihilation velocity setting")
        return capability
    except (OSError, subprocess.SubprocessError, ValueError, KeyError, TypeError, AttributeError) as error:
        raise ValueError(
            f"Planck CMB requires a rebuilt TRSM driver supporting --capabilities "
            f"and --planck-cmb: {path} ({error})"
        ) from error


def require_cmb_capability(executable):
    path = Path(executable).expanduser().resolve()
    stat = path.stat()
    # Return a copy; cached capability data must not be mutated by callers.
    return dict(_driver_capability(str(path), stat.st_mtime_ns, stat.st_size))


@dataclass(frozen=True)
class CMBSignal:
    available: bool = False
    ratio_raw: float = math.nan
    status: str = "missing_output"
    reason: str = "The driver did not report a CMB result"


def parse_cmb_signal(text):
    lines = [line[len(CMB_OUTPUT_PREFIX):] for line in text.splitlines() if line.startswith(CMB_OUTPUT_PREFIX)]
    if not lines:
        return CMBSignal()
    try:
        if len(lines) != 1:
            raise ValueError("duplicate CMB result")
        data = json.loads(lines[0])
        if data["status"] != "ok":
            return CMBSignal(status="calculation_error", reason=str(data.get("reason", "CMB calculation failed")))
        ratio = data["ratio_raw"]
        rate = data["sigma_v_cm3_s"]
        for value in (ratio, rate):
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError("CMB ratio and annihilation rate must be finite and non-negative")
        return CMBSignal(True, float(ratio), "ok", "")
    except (ValueError, KeyError, TypeError, AttributeError) as error:
        return CMBSignal(status="invalid_output", reason=str(error))


@dataclass(frozen=True)
class CMBLimitResult:
    signal: CMBSignal
    fraction: float
    ratio: float
    excluded: bool | None

    @property
    def passed(self):
        return self.signal.available and self.excluded is False


def assess_cmb_limit(signal, omega, rescale=True):
    if not math.isfinite(omega) or omega < 0:
        raise ValueError("CMB rescaling requires a finite non-negative relic density")
    fraction = min(1.0, omega / CMB_OMEGA_REFERENCE) if rescale else 1.0
    if signal.available and (not math.isfinite(signal.ratio_raw) or signal.ratio_raw < 0):
        signal = CMBSignal(status="invalid_output", reason="CMB ratio must be finite and non-negative")
    if not signal.available:
        return CMBLimitResult(signal, fraction, math.nan, None)
    ratio = signal.ratio_raw * fraction**2
    return CMBLimitResult(signal, fraction, ratio, ratio > 1.0)


def cmb_diagnostics(result=None, *, enabled=False, reason="DM calculation unavailable"):
    if result is None:
        return dict(zip(CMB_COLUMNS, (
            enabled, False, "dm_unavailable" if enabled else "disabled",
            reason if enabled else "", None, None, None, None,
        )))
    return dict(zip(CMB_COLUMNS, (
        True, result.signal.available, result.signal.status, result.signal.reason,
        result.signal.ratio_raw if result.signal.available else None,
        result.fraction, result.ratio if result.signal.available else None,
        result.excluded,
    )))
