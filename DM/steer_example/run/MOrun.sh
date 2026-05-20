#!/usr/bin/env bash

shopt -s nullglob

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

MICROMEGAS_MAIN=${MICROMEGAS_MAIN:-$SCRIPT_DIR/../../micromegas_6.1.15/TRSM/main}
EXCLUDER=${EXCLUDER:-"$SCRIPT_DIR/../source/mO_excluder"}
CARD_DIR=${CARD_DIR:-"$SCRIPT_DIR/cards"}
OUTPUT_DIR=${OUTPUT_DIR:-"$SCRIPT_DIR/../output"}
OKS_FILE=${OKS_FILE:-"$SCRIPT_DIR/oks.dat"}

extract_input_value() {
  awk -v key="$1" '$1 == key { print $2; exit }' "$2"
}

extract_omega() {
  awk '
    match($0, /Omega=([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)/, result) {
      print result[1]
      exit
    }
  ' "$1"
}

extract_neutron_si_cross_section() {
  awk '
    /^[[:space:]]*neutron[[:space:]]+SI[[:space:]]/ {
      print $3
      exit
    }
  ' "$1"
}

extract_fermi_lat_line_args() {
  awk '
    function parse_value(line, key, value) {
      if (match(line, key "[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?", value)) {
        sub(key, "", value[0])
        return value[0]
      }
      return ""
    }

    /FermiLAT_line_channel/ {
      e_gamma = parse_value($0, "E_gamma=")
      phi_r16 = parse_value($0, "Phi_R16=")
      if (e_gamma != "" && phi_r16 != "") {
        print e_gamma
        print phi_r16
      }
    }
  ' "$1"
}

safe_move() {
  if [ -e "$1" ]; then
    mv "$1" "$2"
  fi
}

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/OUT_mO"
mkdir -p "$OUTPUT_DIR/DM_data"
mkdir -p "$OUTPUT_DIR/DM_EXCLUDED"
mkdir -p "$OUTPUT_DIR/RelDens_EXCLUDED"
mkdir -p "$OUTPUT_DIR/DirDet_EXCLUDED"
mkdir -p "$OUTPUT_DIR/IndirDet_EXCLUDED"

export MO_EXCLUDER_OUTPUT_DIR="$OUTPUT_DIR"
export MO_EXCLUDER_OKS_FILE="$OKS_FILE"

LOG_FILE=${LOG_FILE:-"$OUTPUT_DIR/MOrun.log"}
if [ -z "${MORUN_LOGGING_INITIALIZED:-}" ]; then
  export MORUN_LOGGING_INITIALIZED=1
  export LOG_FILE
  exec "${BASH:-/bin/bash}" "$0" "$@" 2>&1 | tee -a "$LOG_FILE"
  exit $?
fi

input_files=("$CARD_DIR"/MO_inp*.dat)

if [ ${#input_files[@]} -eq 0 ]; then
  echo "No MO_inp*.dat files found in $CARD_DIR"
  exit 0
fi

if [ ! -x "$MICROMEGAS_MAIN" ]; then
  echo "micrOMEGAs executable not found: $MICROMEGAS_MAIN" >&2
  exit 1
fi

for input_file in "${input_files[@]}"; do
  input_name=${input_file##*/}
  index=${input_name#MO_inp}
  index=${index%.dat}

  echo "$index $input_file"
  echo "Processing $input_file with micromegas"

  "$MICROMEGAS_MAIN" "$input_file" > MO_out

  # MDM=$(extract_mdm MO_out) 
  MDM=$(extract_input_value "MX" "$input_file")
  omega=$(extract_omega MO_out)
  dirdet=$(extract_neutron_si_cross_section MO_out)
  mapfile -t fermi_lat_line_args < <(extract_fermi_lat_line_args MO_out)
  

  lx=$(extract_input_value "LX" "$input_file")
  lhx=$(extract_input_value "LHX" "$input_file")
  lsx=$(extract_input_value "LSX" "$input_file")
  mx=$(extract_input_value "MX" "$input_file")
  vevs=$(extract_input_value "vevs" "$input_file")
  sint=$(extract_input_value "SinT" "$input_file")
  mh2=$(extract_input_value "Mh2" "$input_file")

  echo "$MDM $omega $dirdet"
  echo "$index $lx $lhx $lsx $mx $vevs $sint $mh2 $MDM $omega $dirdet" >> "$OUTPUT_DIR/scan_results.dat"

  if [ -x "$EXCLUDER" ] && [ -f "$OKS_FILE" ]; then
    (
      cd "$SCRIPT_DIR/../source" || exit 1
      MO_EXCLUDER_OUTPUT_DIR="$OUTPUT_DIR" MO_EXCLUDER_OKS_FILE="$OKS_FILE" "$EXCLUDER" "$index" "$MDM" "$omega" "$dirdet" "${fermi_lat_line_args[@]}"
    )
  else
    echo "Skipping mO_excluder for $input_file: $EXCLUDER or $OKS_FILE missing"
  fi
  
  safe_move MO_out "$OUTPUT_DIR/OUT_mO/OUT_mO_${index}"
  safe_move "$OUTPUT_DIR/DM_data/DM_data" "$OUTPUT_DIR/DM_data/DM_data_${index}"
  safe_move "$OUTPUT_DIR/DM_EXCLUDED/DM_EXCLUDED" "$OUTPUT_DIR/DM_EXCLUDED/DM_EXCLUDED_${index}"
  safe_move "$OUTPUT_DIR/RelDens_EXCLUDED/RelDens_EXCLUDED" "$OUTPUT_DIR/RelDens_EXCLUDED/RelDens_EXCLUDED_${index}"
  safe_move "$OUTPUT_DIR/DirDet_EXCLUDED/DirDet_EXCLUDED" "$OUTPUT_DIR/DirDet_EXCLUDED/DirDet_EXCLUDED_${index}"
  safe_move "$OUTPUT_DIR/IndirDet_EXCLUDED/IndirDet_EXCLUDED" "$OUTPUT_DIR/IndirDet_EXCLUDED/IndirDet_EXCLUDED_${index}"

  if [ -f "$OUTPUT_DIR/RelDens_EXCLUDED/RelDens_EXCLUDED_${index}" ]; then
    echo "$index $lx $lhx $lsx $mx $vevs $sint $mh2 $MDM $omega $dirdet" >> "$OUTPUT_DIR/RelDensEXCL.dat"
  fi

  if [ -f "$OUTPUT_DIR/DirDet_EXCLUDED/DirDet_EXCLUDED_${index}" ]; then
    echo "$index $lx $lhx $lsx $mx $vevs $sint $mh2 $MDM $omega $dirdet" >> "$OUTPUT_DIR/DirDetEXCL.dat"
  fi

  if [ -f "$OUTPUT_DIR/IndirDet_EXCLUDED/IndirDet_EXCLUDED_${index}" ]; then
    echo "$index $lx $lhx $lsx $mx $vevs $sint $mh2 $MDM $omega $dirdet" >> "$OUTPUT_DIR/IndirDetEXCL.dat"
  fi
done
