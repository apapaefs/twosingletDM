"""Version selection for the side-by-side TRSM micrOMEGAs installations."""

from pathlib import Path
import os


DEFAULT_MICROMEGAS_VERSION = "7.1.4"
MICROMEGAS_VERSIONS = ("6.1.15", "7.1.4")


def normalize_micromegas_version(value):
    aliases = {"6": "6.1.15", "7": "7.1.4"}
    version = aliases.get(value, value)
    if version not in MICROMEGAS_VERSIONS:
        raise ValueError(f"Unsupported micrOMEGAs version: {value!r}; use 6 or 7")
    return version


def default_micromegas_main(version=DEFAULT_MICROMEGAS_VERSION):
    version = normalize_micromegas_version(version)
    root=Path(os.environ.get("TRSM_RUNTIME_ROOT",Path(__file__).resolve().parents[1] / "runtime-v2"))
    return root / f"micromegas_{version}" / "TRSM" / "main"


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


from functools import lru_cache
import json
import subprocess

@lru_cache(maxsize=16)
def _capabilities(path, mtime, size, virtual_off):
    result=subprocess.run([path,'--capabilities'],check=True,capture_output=True,text=True,timeout=10)
    try: payload=json.loads(result.stdout)
    except ValueError as error: raise ValueError('A rebuilt TRSM driver with v2 capabilities is required: '+path) from error
    if payload.get('physics_version')!='trsm_constraints_v2' or payload.get('loop_hook') is not True:
        raise ValueError('Constraint v2 requires a rebuilt TRSM driver with the corrected loop hook: '+path)
    return payload


def require_v2_capability(executable):
    path=Path(executable).resolve();stat=path.stat()
    return _capabilities(str(path),stat.st_mtime_ns,stat.st_size,"TRSM_LEGACY_VIRTUAL_OFF" in os.environ)
