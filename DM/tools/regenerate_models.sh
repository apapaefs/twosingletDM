#!/bin/sh
# Generate both archived CalcHEP variants from the single LanHEP source.
set -eu
lhep=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
root=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
python3 "$root/DM/tools/sync_sm_inputs.py"
temporary=$(mktemp -d "${TMPDIR:-/tmp}/trsm-lanhep.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM
for mode in On Off; do
    mkdir "$temporary/$mode"
    (cd "$temporary/$mode" && "$lhep" -ca -key "h4G=$mode" "$root/DM/models/lanhep_mdl/TRSM_mixed.mdl") > "$temporary/$mode/generation.log" 2>&1
    cat "$temporary/$mode/generation.log"
    if grep -Eiq '(^Error|^Fatal|Error in|error:)' "$temporary/$mode/generation.log"; then
        exit 1
    fi
    test -f "$temporary/$mode/lgrng1.mdl"
    cp "$temporary/$mode/"*1.mdl "$root/DM/models/h4G$mode/"
done
