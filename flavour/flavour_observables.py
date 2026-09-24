from math import sqrt, pi

if __package__:
    from .flavour_inputs import PHYSICS_INPUTS, central
else:  # Retain the collaborator's standalone example.
    from flavour_inputs import PHYSICS_INPUTS, central
from functools import lru_cache
from pathlib import Path
import numpy as np



def FbS(m_hi: float) -> float:
    """
    QCD and bottomonium correction factor for

        Upsilon(1S) -> gamma h_i.

    Parameters
    ----------
    m_hi:
        Scalar mass in GeV.

    Returns
    -------
    float
        Dimensionless correction factor.
    """

    m_upsilon = central(
        PHYSICS_INPUTS["bottomonium"]["1S"]["mass"]
    )

    if m_hi < 0.0:
        raise ValueError("m_hi must be non-negative.")

    if m_hi >= m_upsilon:
        return 0.0

    return (2.0 / 3.0) * (
        1.0 - (m_hi / m_upsilon)**6
    )


def BRUpsgammahiNorm(
    m_hi: float,
    R_ih: float | complex,
) -> float:
    """
    Calculate

        BR(Upsilon(1S) -> gamma h_i)
        --------------------------------.
        BR(Upsilon(1S) -> e+ e-)

    Parameters
    ----------
    m_hi:
        Scalar mass in GeV.

    R_ih:
        Higgs-doublet component of the scalar h_i.
        It may be real or complex.

    Returns
    -------
    float
        Dimensionless normalized branching ratio.
        Returns zero if on-shell production is
        kinematically forbidden.
    """

    GF = central(
        PHYSICS_INPUTS["constants"]["GF"]
    )

    alpha_em = central(
        PHYSICS_INPUTS["constants"]["alpha_em"]
    )

    m_b = central(
        PHYSICS_INPUTS["constants"]["m_b"]
    )

    m_upsilon = central(
        PHYSICS_INPUTS["bottomonium"]["1S"]["mass"]
    )

    if m_hi < 0.0:
        raise ValueError("m_hi must be non-negative.")

    # On-shell Upsilon(1S) -> gamma h_i is forbidden here.
    if m_hi >= m_upsilon:
        return 0.0

    mixing_squared = abs(R_ih)**2

    phase_space = 1.0 - (m_hi / m_upsilon)**2

    correction = FbS(m_hi)

    return (
        GF * m_b**2
        / (sqrt(2.0) * pi * alpha_em)
        * mixing_squared
        * phase_space
        * correction
    )



def BRUpsgammahi(
    m_hi: float,
    R_ih: float | complex,
) -> float:
    """
    Calculate BR(Upsilon(1S) -> gamma h_i).
    """

    BR_upsilon_ee = central(
        PHYSICS_INPUTS["bottomonium"]["1S"]["BR_ee"]
    )

    return (
        BR_upsilon_ee
        * BRUpsgammahiNorm(
            m_hi=m_hi,
            R_ih=R_ih,
        )
    )

################################################################################
################################################################################
################################################################################

# Directory containing the digitized experimental curves.
LIMIT_DATA_DIR = Path(__file__).resolve().parent
MUMU_MIN = 0.212
TAUTAU_MIN = 3.55
SEARCH_MAX = 9.2
LIMIT_FILES = (
    "mumu_B1B3_Belle_digitized.csv", "mumu_B1B3_BaBar_digitized.csv",
    "tau_tau_B1B2_Belle_digitized.csv", "tau_tau_B1B2_BaBar_digitized.csv",
)


@lru_cache(maxsize=None)
def _load_limit_curve(filename: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Load a digitized mass-dependent upper-limit curve.

    The CSV must contain:
        column 1: scalar mass in GeV
        column 2: upper limit on the product branching fraction
    """
    path = LIMIT_DATA_DIR / filename

    data = np.genfromtxt(
        path,
        delimiter=",",
        comments="#",
        skip_header=1,
        dtype=float,
    )

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 2:
        raise ValueError(
            f"{path} must contain at least two columns: "
            "mass_GeV and upper_limit."
        )

    masses = data[:, 0]
    limits = data[:, 1]

    valid = (
        np.isfinite(masses)
        & np.isfinite(limits)
        & (limits > 0.0)
    )

    masses = masses[valid]
    limits = limits[valid]

    if masses.size < 2:
        raise ValueError(f"{path} does not contain enough valid points.")

    # Ensure increasing mass order.
    order = np.argsort(masses)
    masses = masses[order]
    limits = limits[order]

    return masses, limits


def _interpolate_limit(
    m_hi: float,
    filename: str,
) -> float:
    """
    Return the mass-dependent upper limit.

    Interpolation is linear in mass and logarithmic in the upper limit.
    Outside the experimentally covered mass range, return infinity,
    corresponding to no constraint from that experiment.
    """
    if not np.isfinite(m_hi):
        raise ValueError("m_hi must be finite.")

    masses, limits = _load_limit_curve(filename)

    if m_hi < masses[0] or m_hi > masses[-1]:
        return float("inf")

    log_limit = np.interp(
        m_hi,
        masses,
        np.log10(limits),
    )

    return float(10.0**log_limit)


def belle_mumu_limit(m_hi: float) -> float:
    """
    Belle 90% upper limit on

        BR(Upsilon(1S) -> gamma h_i)
        * BR(h_i -> mu+ mu-).
    """
    return _interpolate_limit(
        m_hi,
        "mumu_B1B3_Belle_digitized.csv",
    )


def babar_mumu_limit(m_hi: float) -> float:
    """
    BaBar 90% upper limit on

        BR(Upsilon(1S) -> gamma h_i)
        * BR(h_i -> mu+ mu-).
    """
    return _interpolate_limit(
        m_hi,
        "mumu_B1B3_BaBar_digitized.csv",
    )


def belle_tautau_limit(m_hi: float) -> float:
    """
    Belle 90% upper limit on

        BR(Upsilon(1S) -> gamma h_i)
        * BR(h_i -> tau+ tau-).

    Returns infinity outside the mass range covered by the CSV data.
    """
    return _interpolate_limit(
        m_hi=m_hi,
        filename="tau_tau_B1B2_Belle_digitized.csv",
    )


def babar_tautau_limit(m_hi: float) -> float:
    """
    BaBar 90% upper limit on

        BR(Upsilon(1S) -> gamma h_i)
        * BR(h_i -> tau+ tau-).

    Returns infinity outside the mass range covered by the CSV data.
    """
    return _interpolate_limit(
        m_hi=m_hi,
        filename="tau_tau_B1B2_BaBar_digitized.csv",
    )




##################################################################################
##################################################################################
##################################################################################


def BRUpsgammahiFinalState(
    m_hi: float,
    R_ih: float | complex,
    BR_hi_to_final: float,
) -> float:
    """
    Calculate

        BR(Upsilon(1S) -> gamma h_i)
        x BR(h_i -> final state).

    This is the quantity constrained by Belle and BaBar.

    Parameters
    ----------
    m_hi:
        Scalar mass in GeV.

    R_ih:
        Higgs-doublet component of h_i.

    BR_hi_to_final:
        Branching fraction of h_i into the selected final state,
        for example h_i -> mu+ mu- or h_i -> tau+ tau-.

    Returns
    -------
    float
        Product branching fraction.
    """

    if not 0.0 <= BR_hi_to_final <= 1.0:
        raise ValueError(
            "BR_hi_to_final must lie between zero and one."
        )

    production_BR = BRUpsgammahi(
        m_hi=m_hi,
        R_ih=R_ih,
    )

    return production_BR * BR_hi_to_final


def CheckUpsLeptonBounds(
    m_hi: float,
    R_ih: float | complex,
    BR_hi_to_mumu: float,
    BR_hi_to_tautau: float,
) -> int:
    """
    Check Upsilon(1S) radiative-decay constraints for

        h_i -> mu+ mu-
        h_i -> tau+ tau-.

    Returns
    -------
    int
        1 if all applicable constraints are satisfied.
        0 if at least one applicable constraint is violated.
    """

    return int(EvaluateUpsLeptonBounds(
        m_hi, R_ih, BR_hi_to_mumu, BR_hi_to_tautau,
    )["passed"])


def EvaluateUpsLeptonBounds(m_hi, R_ih, BR_hi_to_mumu, BR_hi_to_tautau):
    """Detailed version of the same test; assess both channels without short circuiting.

    Infinite limits are represented by None in diagnostics for portable JSON.
    The stronger available individual limit is used, not a combined likelihood.
    """
    if not np.isfinite(m_hi) or m_hi <= 0.0:
        raise ValueError("m_hi must be finite and positive.")

    if not np.isfinite(abs(R_ih)):
        raise ValueError("R_ih must be finite.")

    if (
        not np.isfinite(BR_hi_to_mumu)
        or not 0.0 <= BR_hi_to_mumu <= 1.0
    ):
        raise ValueError(
            "BR_hi_to_mumu must be finite and lie between 0 and 1."
        )

    if (
        not np.isfinite(BR_hi_to_tautau)
        or not 0.0 <= BR_hi_to_tautau <= 1.0
    ):
        raise ValueError(
            "BR_hi_to_tautau must be finite and lie between 0 and 1."
        )

    production = BRUpsgammahi(m_hi, R_ih)
    channels = {}
    for name, minimum, br, belle, babar in (
        ("mumu", MUMU_MIN, BR_hi_to_mumu, belle_mumu_limit, babar_mumu_limit),
        ("tautau", TAUTAU_MIN, BR_hi_to_tautau, belle_tautau_limit, babar_tautau_limit),
    ):
        in_search = minimum <= m_hi <= SEARCH_MAX
        limits = {"Belle": belle(m_hi), "BaBar": babar(m_hi)} if in_search else {}
        available = {k: float(v) for k, v in limits.items() if np.isfinite(v)}
        strongest = min(available, key=available.get) if available else None
        limit = available.get(strongest)
        prediction = production * br
        ratio = float(prediction / limit) if limit is not None else None
        excluded = limit is not None and prediction > limit
        channels[name] = {
            "br": float(br), "prediction": float(prediction),
            "in_search_range": in_search, "covered": bool(available),
            "limits": {key: available.get(key) for key in ("Belle", "BaBar")},
            "limit": limit, "experiment": strongest, "ratio": ratio,
            "excluded": bool(excluded),
            "status": ("outside_search" if not in_search else "outside_curve" if not available
                       else "excluded" if excluded else "zero_signal" if prediction == 0 else "passed"),
        }
    return {"passed": not any(c["excluded"] for c in channels.values()),
            "production_br": float(production), "channels": channels}
