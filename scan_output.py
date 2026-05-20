from pathlib import Path


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
    "dm_limit_model",
    "dm_rescale",
]


def mg5_columns(mg5xsecs):
    return [f"mg5_xsec_{process}_pb" for process in mg5xsecs.keys()]


def format_output_value(value):
    if value is None:
        return "nan"
    return str(value)


def output_columns(mg5xsecs):
    if mg5xsecs is None:
        mg5xsecs = {}
    return POINT_COLUMNS + DM_EXCLUSION_COLUMNS + mg5_columns(mg5xsecs)


def output_row(point_info, mg5xsecs=None):
    if mg5xsecs is None:
        mg5xsecs = {}

    values = [point_info.get(column) for column in POINT_COLUMNS + DM_EXCLUSION_COLUMNS]
    values.extend(mg5xsecs[process] for process in mg5xsecs.keys())
    return "\t".join(format_output_value(value) for value in values)


def write_valid_point(outfile, point_info, mg5xsecs=None):
    outfile = Path(outfile)
    outfile.parent.mkdir(parents=True, exist_ok=True)
    write_header = not outfile.exists() or outfile.stat().st_size == 0

    with outfile.open("a", encoding="ascii") as filestream:
        if write_header:
            filestream.write("\t".join(output_columns(mg5xsecs)) + "\n")
        filestream.write(output_row(point_info, mg5xsecs) + "\n")
