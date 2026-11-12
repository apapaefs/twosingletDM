"""Stored-assessment selection and factorized associated-production rates.

The legacy mono-Z column remains the tree contribution.  The total is the
tree result plus the separate loop-induced gg result, not a full NLO result.
"""

import math

ASSOCIATED_PROCESSES = ("gg_heta0", "pp_eta0Z", "gg_eta0Z")
ADDITIONAL_RATE_COLUMNS = (
    "mg5_xsec_pp_eta0Z_total_pb", "mono_z_gg_xsec_pb", "mono_z_total_xsec_pb",
)
DERIVED_RATE_COLUMNS = ("mono_higgs_xsec_pb", "mono_z_xsec_pb", *ADDITIONAL_RATE_COLUMNS)


def finite_number(value):
    try:
        number = float(value)
    except (ValueError, TypeError):
        return None
    return number if math.isfinite(number) else None


def stored_full_viability(row):
    """Apply the plot suite's version-aware full selection to stored verdicts."""
    def passed(key):
        value = row.get(key)
        return value is True or (isinstance(value, str) and value.strip() == "True")

    if row.get("constraint_version") == "trsm_constraints_v2":
        non_dm = passed("thc") and passed("experimental_subset")
    else:
        non_dm = all(passed(key) for key in ("evo", "thc", "hb", "hs", "ewpo", "wmass"))
    return non_dm and passed("flavour") and passed("dm")


def derive_mg5_rates(values):
    """Return computable rates; unavailable components never count as zero."""
    production = {}
    for process in ASSOCIATED_PROCESSES:
        value = finite_number(values.get(f"mg5_xsec_{process}_pb"))
        if value is not None and value >= 0:
            production[process] = value
    result = {}
    if "pp_eta0Z" in production and "gg_eta0Z" in production:
        total = production["pp_eta0Z"] + production["gg_eta0Z"]
        if math.isfinite(total):
            result["mg5_xsec_pp_eta0Z_total_pb"] = total
    br = finite_number(values.get("h2_h3h3_br"))
    if br is None or not 0 <= br <= 1:
        return result
    for process, column in (("gg_heta0", "mono_higgs_xsec_pb"),
                            ("pp_eta0Z", "mono_z_xsec_pb"),
                            ("gg_eta0Z", "mono_z_gg_xsec_pb")):
        if process in production:
            result[column] = production[process] * br
    if "mg5_xsec_pp_eta0Z_total_pb" in result:
        result["mono_z_total_xsec_pb"] = result["mg5_xsec_pp_eta0Z_total_pb"] * br
    return result
