"""Published, elastic SI upper-limit tables used after the micrOMEGAs run."""

import hashlib
import json
import math
import re
from bisect import bisect_left
from dataclasses import dataclass
from pathlib import Path


DEFAULT_LIMIT_MODEL = "lz2025-source"
SI_TABLE_SCHEMA = "trsm_si_upper_limit_v1"


@dataclass(frozen=True)
class SILimitTable:
    label: str
    source: str
    path: Path
    sha256: str
    input_unit: str
    masses_gev: tuple[float, ...]
    limits_pb: tuple[float, ...]
    provenance_json: str | None = None

    @property
    def model_id(self):
        return f"table:{self.label}:{self.sha256[:12]}"

    def upper_limit_pb(self, mass_gev):
        """Interpolate in log(mass), log(cross section); never extrapolate."""
        if not math.isfinite(mass_gev) or mass_gev <= 0:
            raise ValueError("Dark matter mass must be finite and positive")
        if not self.masses_gev[0] <= mass_gev <= self.masses_gev[-1]:
            raise ValueError(
                f"DM mass {mass_gev:g} GeV is outside the {self.label} table "
                f"range [{self.masses_gev[0]:g}, {self.masses_gev[-1]:g}] GeV; "
                "extrapolation is disabled"
            )
        index = bisect_left(self.masses_gev, mass_gev)
        if self.masses_gev[index] == mass_gev:
            return self.limits_pb[index]
        low, high = self.masses_gev[index - 1:index + 1]
        weight = math.log(mass_gev / low) / math.log(high / low)
        return math.exp(
            (1 - weight) * math.log(self.limits_pb[index - 1])
            + weight * math.log(self.limits_pb[index])
        )

    def metadata(self):
        metadata = {
            "model": self.model_id,
            "table_path": str(self.path),
            "sha256": self.sha256,
            "source": self.source,
            "confidence_level": 0.9,
            "interaction": "elastic_isoscalar_si",
            "limit_kind": "observed_upper",
            "cross_section": "per_nucleon",
            "input_unit": self.input_unit,
            "comparison_unit": "pb",
            "mass_range_gev": [self.masses_gev[0], self.masses_gev[-1]],
            "interpolation": "log-log",
            "extrapolation": "error",
        }
        if self.provenance_json is not None:
            metadata["provenance"] = json.loads(self.provenance_json)
        return metadata


def _positive_number(value, label):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise ValueError(f"{label} must be a finite, positive number")
    return float(value)


def load_si_limit_table(path):
    """Load a documented 90% observed upper limit, with explicit physical units.

    Input is the normalized JSON format documented in DM/direct-detection.md,
    not an arbitrary HEPData export or a table of EFT coupling coefficients.
    The immutable object freezes both numbers and provenance for the scan.
    """
    path = Path(path).expanduser().resolve()
    raw = path.read_bytes()
    data = json.loads(raw)
    if not isinstance(data, dict) or data.get("schema") != SI_TABLE_SCHEMA:
        raise ValueError(f"SI limit table must use schema {SI_TABLE_SCHEMA}")
    for key, expected in (
        ("confidence_level", 0.9),
        ("interaction", "elastic_isoscalar_si"),
        ("limit_kind", "observed_upper"),
        ("cross_section", "per_nucleon"),
    ):
        if data.get(key) != expected:
            raise ValueError(f"SI limit table requires {key}={expected!r}")
    label = data.get("label", "")
    if not isinstance(label, str) or not re.fullmatch(r"[a-z0-9][a-z0-9._-]{0,63}", label):
        raise ValueError("SI table label must contain 1-64 lowercase letters, digits, '.', '_' or '-'")
    source = data.get("source")
    if not isinstance(source, str) or not source.strip():
        raise ValueError("SI table requires a nonempty publication/data source")
    provenance = data.get("provenance")
    if provenance is not None and not isinstance(provenance, dict):
        raise ValueError("SI table provenance must be an object")
    provenance_json = (
        json.dumps(provenance, sort_keys=True, allow_nan=False)
        if provenance is not None else None
    )
    unit = data.get("cross_section_unit")
    if unit not in ("pb", "cm2"):
        raise ValueError("SI table cross_section_unit must be 'pb' or 'cm2'")
    points = data.get("points")
    if not isinstance(points, list) or len(points) < 2:
        raise ValueError("SI table requires at least two mass/upper-limit points")
    masses, limits = [], []
    for index, point in enumerate(points):
        if not isinstance(point, dict):
            raise ValueError(f"SI table point {index} must be an object")
        mass = _positive_number(point.get("mass_GeV"), f"Point {index} mass_GeV")
        limit = _positive_number(point.get("upper_limit"), f"Point {index} upper_limit")
        if masses and mass <= masses[-1]:
            raise ValueError("SI table masses must be strictly increasing")
        masses.append(mass)
        # 1 pb = 10^-36 cm^2. micrOMEGAs prints nucleon cross sections in pb.
        limits.append(_positive_number(limit * (1e36 if unit == "cm2" else 1), "Limit in pb"))
    return SILimitTable(
        label, source.strip(), path, hashlib.sha256(raw).hexdigest(), unit,
        tuple(masses), tuple(limits), provenance_json,
    )


def direct_detection_configuration(args):
    table = getattr(args, "_dm_limit_table", None)
    if table is None:
        return {"model": DEFAULT_LIMIT_MODEL}
    return table.metadata()
