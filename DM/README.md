# Dark matter: installation and standalone scans

Start with a full clone of `twosingletDM`. For a DM-only scan you need Python
3.10+, a C/C++ compiler, a Fortran compiler supported by micrOMEGAs, `make`,
`patch` and `tar`. NumPy and matplotlib are needed for plotting. HiggsTools,
its datasets and BSMPT are not needed for `DM/steer_example`.

## Install micrOMEGAs

Download the **7.1.4** source archive from the
[micrOMEGAs distribution](https://micromegasdm.github.io/). Keep older working
installations intact. From the repository root:

```bash
# Set this to the archive you downloaded.
TRSM_MO_ARCHIVE="$HOME/Downloads/micromegas_7.1.4.tgz"
mkdir -p ../runtime-v2
tar -xzf "$TRSM_MO_ARCHIVE" -C ../runtime-v2
sh DM/setup_micromegas.sh ../runtime-v2/micromegas_7.1.4

python3 -m venv .venv-dm
source .venv-dm/bin/activate
python -m pip install numpy matplotlib
python DM/steer_example/check_installation.py
```

Extract only into a location where `micromegas_7.1.4` does not already exist.
The setup script patches the spectrum finite-value check, builds micrOMEGAs
serially, creates `TRSM`, and installs `main.c`, `trsm_loop.c`, the generated
`models/h4GOn/*.mdl` files and `data.par`, then compiles `TRSM/main`. It refuses
to replace an existing `TRSM` installation. The released archive used for v2
has SHA256 `c8cf207b17541a5b7d7e7ff157f1c1eb36e1dd3f7547e7745559b195c2a34264`.

For a different destination, extract there and pass the resulting release
directory to setup. Then supply `--micromegas-main /absolute/path/TRSM/main`
to the Python runner/checker, or export `MICROMEGAS_MAIN`.
`TRSM_RUNTIME_ROOT` selects a common parent for side-by-side installations.
The same procedure accepts the **6.1.15** archive/directory; run with
`--micromegas-version 6` to select that backend.

`source setup_micromegas.sh` without an argument is not an installation command:
setup requires the extracted directory. It now prints help and returns status 2
without closing the caller or changing its shell options/directory. Prefer
`sh setup_micromegas.sh /path/to/micromegas_7.1.4` when already inside `DM`.
`--help` prints usage without building anything.

If compilation fails, retain the first compiler error from the build log.
Consult the extracted release's `README` for compiler selection. The setup
script does not download compilers or change shell profiles. Use a fresh
extraction after a failed partial installation instead of overwriting a
possibly active model.

## Run and inspect

See [the general scan, card, shell and plotting walkthrough](steer_example/README.md).
It starts from `DM/steer_example`, preserves the old input and output layouts,
and explains the optional CMB checks separately. A compiled capability check is:

```bash
../runtime-v2/micromegas_7.1.4/TRSM/main --capabilities
```

The Python runners create a private writable CalcHEP directory automatically.
For a direct diagnostic invocation of the native executable, set a fresh
`TRSM_RUNTIME_DIR` yourself; this is unnecessary through the Python/shell scan
entry points.

- [What changed in the model, including `x1`](models/changes-v2.md)
- [Planck CMB prescription](planck-cmb.md)
- [Direct-detection table and conventions](direct-detection.md)
- [Full constraint profile and runtime setup](../docs/constraints-v2.md)

The default backend is micrOMEGAs 7.1.4 with the v2 model and loop helper.
Changing interfaces back to their historical form does not restore older
physics defaults: the model, limits and abundance conventions are recorded
explicitly in each run's metadata.
