#!/bin/sh
# Configure a freshly extracted micrOMEGAs release with the canonical TRSM model.
set -eu

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/micromegas_6.1.15-or-7.1.4" >&2
    exit 2
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
install_dir=$(CDPATH= cd -- "$1" && pwd)
version=${install_dir##*/micromegas_}
case "$version" in
    6.1.15|7.1.4) ;;
    *) echo "Unsupported micrOMEGAs installation: $install_dir" >&2; exit 2 ;;
esac

# Never rebuild a model that may be in use by a campaign.
if [ -e "$install_dir/TRSM" ]; then
    echo "Refusing to overwrite existing $install_dir/TRSM; use a fresh extraction." >&2
    exit 2
fi

cd "$install_dir"
patch -p1 < "$script_dir/patches/micromegas-$version-calcspectrum-finite-guard.patch"
# Upstream Makefiles update shared static archives; build serially.
make
./newProject TRSM
cp "$script_dir/main.c" TRSM/main.c
cp "$script_dir/trsm_loop.c" TRSM/lib/trsm_loop.c
if [ "$version" = 7.1.4 ]; then
    sed '1i\
#define TRSM_MO7 1\
' "$script_dir/trsm_loop.c" > TRSM/lib/trsm_loop.c
fi
cp "$script_dir"/models/h4GOn/*.mdl TRSM/work/models/
printf '\nextern double trsm_loop_abs(double,double);\n' >> TRSM/work/models/extlib1.mdl
cp "$script_dir/data.par" TRSM/data.par
make -C TRSM main=main.c
echo "Installed micrOMEGAs $version: $install_dir/TRSM/main"
