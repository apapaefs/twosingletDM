#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
args=(--cards "${CARD_DIR:-$SCRIPT_DIR/cards}" --output-dir "${OUTPUT_DIR:-$SCRIPT_DIR/../output}")
if [[ -n "${MICROMEGAS_MAIN:-}" ]]; then
    args+=(--micromegas-main "$MICROMEGAS_MAIN")
fi
if [[ -n "${LOG_FILE:-}" ]]; then
    args+=(--log-file "$LOG_FILE")
fi
exec "${PYTHON:-python3}" "$SCRIPT_DIR/../run_scan.py" "${args[@]}" "$@"
