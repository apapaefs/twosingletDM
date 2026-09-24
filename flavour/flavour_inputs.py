from typing import Any


# All masses and widths are expressed in GeV.
PHYSICS_INPUTS: dict[str, Any] = {
    "metadata": {
        "source": "PDG Live",
        "accessed": "2026-09-24",
    },

    "constants": {
        "GF": {
            "value": 1.1663787e-5,
            "uncertainty": None,
            "unit": "GeV^-2",
        },
        "alpha_em": {
            "value": 1.0 / 137.035999084,
            "uncertainty": None,
            "unit": "dimensionless",
        },
        "m_b": {
            "value": 4.18,
            "uncertainty": 0.03,
            "unit": "GeV",
            "scheme": "MSbar",
            "scale": "m_b",
        },
    },

    "bottomonium": {
        "1S": {
            "mass": {
                "value": 9.46040,
                "uncertainty": 0.00010,
                "unit": "GeV",
                "source": "PDG Live",
            },
            "width": {
                "value": 54.02e-6,
                "uncertainty": 1.25e-6,
                "unit": "GeV",
                "source": "PDG Live",
            },
            "BR_ee": {
                "value": 0.0239,
                "uncertainty": 0.0008,
                "unit": "dimensionless",
                "source": "PDG Live",
            },
            "BR_mumu": {
                "value": 0.0248,
                "uncertainty": 0.0004,
                "unit": "dimensionless",
                "source": "PDG Live",
            },
            "BR_tautau": {
                "value": 0.0260,
                "uncertainty": 0.0010,
                "unit": "dimensionless",
                "source": "PDG Live",
            },
        },

        "2S": {
            "mass": {
                "value": 10.0234,
                "uncertainty": 0.0005,
                "unit": "GeV",
                "source": "PDG Live",
            },
            "width": {
                "value": 31.98e-6,
                "uncertainty": 2.63e-6,
                "unit": "GeV",
                "source": "PDG Live",
            },
            "BR_ee": {
                "value": 0.0191,
                "uncertainty": 0.0016,
                "unit": "dimensionless",
                "source": "PDG Live",
            },
            "BR_mumu": {
                "value": 0.0193,
                "uncertainty": 0.0017,
                "unit": "dimensionless",
                "source": "PDG Live",
            },
            "BR_tautau": {
                "value": 0.0200,
                "uncertainty": 0.0021,
                "unit": "dimensionless",
                "source": "PDG Live",
            },
        },
    },
}




def central(parameter: Any) -> float:
    """
    Return the central numerical value of a physics input.

    Accepted formats
    ----------------
    1. A direct numerical value:
           9.4604

    2. A dictionary containing:
           {"value": 9.4604, ...}

    3. An uncertainties object with a nominal_value attribute.
    """

    if isinstance(parameter, (int, float)):
        return float(parameter)

    if isinstance(parameter, dict):
        if "value" in parameter:
            return float(parameter["value"])

        if "central" in parameter:
            return float(parameter["central"])

        raise KeyError(
            "Physics-input dictionary must contain "
            "a 'value' or 'central' entry."
        )

    if hasattr(parameter, "nominal_value"):
        return float(parameter.nominal_value)

    raise TypeError(
        f"Cannot extract a central value from {parameter!r}."
    )
