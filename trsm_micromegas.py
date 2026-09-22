"""Version selection for the side-by-side TRSM micrOMEGAs installations."""

from pathlib import Path


DEFAULT_MICROMEGAS_VERSION = "6.1.15"
MICROMEGAS_VERSIONS = ("6.1.15", "7.1.4")


def normalize_micromegas_version(value):
    aliases = {"6": "6.1.15", "7": "7.1.4"}
    version = aliases.get(value, value)
    if version not in MICROMEGAS_VERSIONS:
        raise ValueError(f"Unsupported micrOMEGAs version: {value!r}; use 6 or 7")
    return version


def default_micromegas_main(version=DEFAULT_MICROMEGAS_VERSION):
    version = normalize_micromegas_version(version)
    return Path(__file__).resolve().parents[1] / f"micromegas_{version}" / "TRSM" / "main"


def micromegas_configuration(args):
    version = normalize_micromegas_version(
        getattr(args, "micromegas_version", DEFAULT_MICROMEGAS_VERSION)
    )
    executable = getattr(args, "micromegas_main", None)
    if executable is None:
        executable = default_micromegas_main(version)
    return {
        "version": version,
        "executable": str(Path(executable).expanduser().resolve()),
    }
