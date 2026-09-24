from pathlib import Path

from dm_thermal_relic_diagnostic import RESONANCE_COLUMNS, THERMAL_VEV_COLUMNS, THERMAL_COUPLING_COLUMNS
from ewpt_entry_criterion import EW_ENTRY_COLUMNS
from ewpt_x_history import X_HISTORY_COLUMNS
from trsm_cmb import CMB_COLUMNS
from trsm_theory_diagnostics import THEORY_COLUMNS
from trsm_constraint_profile import PROFILE_COLUMNS
from ewpt_assessment import STATUS_COLUMNS
from trsm_flavour import FLAVOUR_COLUMNS

V2_COLUMNS = list(PROFILE_COLUMNS) + list(THEORY_COLUMNS) + list(STATUS_COLUMNS) + [
    "dm_solver_error", "dm_calculation_status", "dm_assessment_reason",
    "dm_direct_detection_available", "dm_actual_inputs", "dm_h1_width_GeV", "dm_h2_width_GeV",
    "dm_abundance_reference", "dm_loop_hook_relic_calls", "dm_loop_hook_indirect_calls",
    "ewpt_constraint_version",
]


POINT_COLUMNS = [
    "M2",
    "M3",
    "vs",
    "vx",
    "a12",
    "a13",
    "a23",
    "lX",
    "lPhiX",
    "lSX",
    "w1",
    "w2",
    "w3",
    "h1_h3h3_width",
    "h1_h3h3_br",
    "h1_h2h2_width",
    "h1_h2h2_br",
    "h2_h3h3_width",
    "h2_h3h3_br",
    "h2_h1h1_width",
    "h2_h1h1_br",
    "xs136_lo_h1_pb",
    "xs136_lo_h2_pb",
    "xsec_h2_h1h1_one_h1_invisible_pb",
    "xsec_h1_h2h2_one_h2_invisible_pb",
    "mono_higgs_xsec_pb",
    "mono_z_xsec_pb",
    "higgs_invisible_widths_included",
    "portal_convention",
    "micromegas_model_convention",
    "k1",
    "k2",
    "k3",
    "K111",
    "K112",
    "K113",
    "K123",
    "K122",
    "K1111",
    "K1112",
    "K1113",
    "K133",
    "K233",
    "evo",
    "thc",
    "hb",
    "hs",
    "ewpo",
    "wmass",
    "dm",
]


DM_EXCLUSION_COLUMNS = [
    "dm_mdm",
    "dm_omega",
    "dm_relic_upper_limit",
    "dm_dir_det",
    "dm_dir_det_limit",
    "dm_lux_base_limit",
    "dm_relic_excluded",
    "dm_direct_detection_excluded",
    "dm_indirect_available",
    "dm_indirect_channels_seen",
    "dm_indirect_channels_used",
    "dm_indirect_energy",
    "dm_indirect_flux",
    "dm_indirect_limit",
    "dm_indirect_ratio",
    "dm_indirect_detection_excluded",
    "dm_limit_model",
    "dm_rescale",
    "dm_xf",
    "dm_freezeout_temperature_GeV",
    *RESONANCE_COLUMNS,
]


HIGGSTOOLS_COLUMNS = [
    "higgstools_hb_selected_limits",
    "higgstools_hb_top_obs",
    "higgstools_hs_chi2",
    "higgstools_hs_delta_chi2",
    "higgstools_hs_top_chi2",
]


EWPT_COLUMNS = [
    "ewpt_ew_true_over_T",
    "ewpt_ew_jump_over_T",
    "ewpt_global_phase_path",
    "ewpt_has_x_broken",
    "ewpt_ew_step_index",
    "ewpt_status",
    "ewpt_error",
    *X_HISTORY_COLUMNS,
    *THERMAL_VEV_COLUMNS,
    *THERMAL_COUPLING_COLUMNS,
]


def mg5_columns(mg5xsecs):
    return [f"mg5_xsec_{process}_pb" for process in mg5xsecs.keys()]


def format_output_value(value):
    if value is None:
        return "nan"
    return str(value)


def output_columns(mg5xsecs, *, planck_cmb=False):
    if mg5xsecs is None:
        mg5xsecs = {}
    cmb_columns = list(CMB_COLUMNS) if planck_cmb else []
    return POINT_COLUMNS + DM_EXCLUSION_COLUMNS + EWPT_COLUMNS + HIGGSTOOLS_COLUMNS + mg5_columns(mg5xsecs) + cmb_columns + list(EW_ENTRY_COLUMNS) + V2_COLUMNS + list(FLAVOUR_COLUMNS)


def output_row(point_info, mg5xsecs=None, *, planck_cmb=False):
    if mg5xsecs is None:
        mg5xsecs = {}

    values = [point_info.get(column) for column in POINT_COLUMNS + DM_EXCLUSION_COLUMNS + EWPT_COLUMNS + HIGGSTOOLS_COLUMNS]
    values.extend(mg5xsecs[process] for process in mg5xsecs.keys())
    if planck_cmb:
        values.extend(point_info.get(column) for column in CMB_COLUMNS)
    values.extend(point_info.get(column) for column in EW_ENTRY_COLUMNS)
    values.extend(point_info.get(column) for column in V2_COLUMNS)
    values.extend(point_info.get(column) for column in FLAVOUR_COLUMNS)
    return "\t".join(format_output_value(value) for value in values)


def write_valid_point(outfile, point_info, mg5xsecs=None, *, planck_cmb=False):
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    write_header = not outfile.exists() or outfile.stat().st_size == 0

    with outfile.open("a", encoding="ascii") as filestream:
        if write_header:
            filestream.write("\t".join(output_columns(mg5xsecs, planck_cmb=planck_cmb)) + "\n")
        filestream.write(output_row(point_info, mg5xsecs, planck_cmb=planck_cmb) + "\n")
