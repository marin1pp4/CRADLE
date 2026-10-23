#!/usr/bin/env bash
# Grid-scan convergence study for CRADLE++ (polarized beta+ and beta-).
# Place this script in the root of cradle_dev and run it from anywhere.
# Usage: bash gridscan_resolution_sweep.sh [angular|energy|all]
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="$ROOT/build"
EXEC="$BUILD/CRADLE++"
CONFIG_DIR="$ROOT/config"
MODE="${1:-all}"

if [[ "$MODE" != "all" && "$MODE" != "angular" && "$MODE" != "energy" ]]; then
  echo "Usage: bash $0 [angular|energy|all]" >&2
  exit 2
fi
if [[ ! -x "$EXEC" ]]; then
  echo "ERROR: Executable not found: $EXEC" >&2
  exit 2
fi

# Match these filenames to your existing configuration files.
PLUS_CONFIG="$CONFIG_DIR/config_gridScan.txt"
MINUS_CONFIG="$CONFIG_DIR/config_gridScan_tritium.txt"
for config in "$PLUS_CONFIG" "$MINUS_CONFIG"; do
  if [[ ! -f "$config" ]]; then
    echo "ERROR: Configuration file not found: $config" >&2
    exit 2
  fi
  for option in GridScan GridScanNE GridScanNU GridScanNV; do
    if ! grep -Eq "^[[:space:]]*${option}[[:space:]]*=" "$config"; then
      echo "ERROR: Missing '${option}=' in $config" >&2
      exit 2
    fi
  done
done

RUN_ID="$(date +%Y%m%d_%H%M%S)"
OUT="$ROOT/gridscan_convergence_$RUN_ID"
mkdir -p "$OUT/logs" "$OUT/events"
CSV="$OUT/results.csv"
printf 'nucleus,decay_mode,sweep,NE,NU,NV,grid_points,W_min,Ee_keV,u,v,scan_seconds,validation,exit_code,log_file\n' > "$CSV"

# Extended angular convergence study. Change the arrays as needed.
# Each run is a separate process and uses one MC thread to reduce unrelated concurrency issues.
ANGULAR_POINTS=(20 30 50 70 100 150 200 250 300 400 600 800 1000)
ENERGY_POINTS=(10 20 30 50 75 100)
ANGULAR_FIXED_NE=50
ENERGY_FIXED_NUV=100

# Extract first field from the formatted grid-scan report, returning NA if absent.
read_field() {
  local label="$1" file="$2" value
  value="$(awk -F ':' -v label="$label" '
    { key=$1; gsub(/^[[:space:]]+|[[:space:]]+$/, "", key) }
    key == label {
      value=$2;
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", value);
      split(value, parts, /[[:space:]]+/);
      print parts[1];
      exit
    }
  ' "$file")"
  printf '%s' "${value:-NA}"
}

# Use original config as a template; never edit it in place.
make_config() {
  local source="$1" dest="$2" ne="$3" nu="$4" nv="$5"
  sed -E \
    -e 's/^[[:space:]]*GridScan[[:space:]]*=.*/GridScan=true/' \
    -e "s/^[[:space:]]*GridScanNE[[:space:]]*=.*/GridScanNE=${ne}/" \
    -e "s/^[[:space:]]*GridScanNU[[:space:]]*=.*/GridScanNU=${nu}/" \
    -e "s/^[[:space:]]*GridScanNV[[:space:]]*=.*/GridScanNV=${nv}/" \
    "$source" > "$dest"
}

# Temporary file resides within config/ for convenient inspection while running.
TMP_CONFIG="$(mktemp "$CONFIG_DIR/.gridscan_sweep_XXXXXXXX.txt")"
trap 'rm -f "$TMP_CONFIG"' EXIT

run_one() {
  local nucleus="$1" decay_mode="$2" source_config="$3" sweep="$4"
  local ne="$5" nu="$6" nv="$7"
  local tag="${nucleus}_${sweep}_E${ne}_U${nu}_V${nv}"
  local logfile="$OUT/logs/${tag}.log"
  local eventfile="$OUT/events/${tag}.txt"
  local exit_code=0 validation wmin ee u v seconds points
  points=$((ne * nu * nv))

  make_config "$source_config" "$TMP_CONFIG" "$ne" "$nu" "$nv"
  echo "[$nucleus | $sweep] NE=$ne NU=$nu NV=$nv ($points points)"

  # Running from build/ preserves existing relative paths to the custom JSON.
  # Keep the radiative-correction CLI override, as requested.
  if (
    cd "$BUILD"
    "$EXEC" \
      -c "$TMP_CONFIG" \
      BetaDecay --radiativecorrection false \
      nucleus -n "$nucleus" -Z "$([[ "$nucleus" == '11C' ]] && echo 6 || echo 1)" \
        -A "$([[ "$nucleus" == '11C' ]] && echo 11 || echo 3)" -e 0 \
      general -l 10 -t 1 -s 2026999 -o "$eventfile"
  ) > "$logfile" 2>&1; then
    exit_code=0
  else
    exit_code=$?
  fi

  wmin="$(read_field 'W_min' "$logfile")"
  ee="$(read_field 'Electron energy' "$logfile")"
  u="$(read_field 'u' "$logfile")"
  v="$(read_field 'v' "$logfile")"
  seconds="$(read_field 'Grid scan time' "$logfile")"
  validation="$(read_field 'Validation' "$logfile")"

  if [[ "$validation" == 'NA' ]]; then
    validation='ERROR_NO_SCAN'
  fi
  printf '%s,%s,%s,%d,%d,%d,%d,%s,%s,%s,%s,%s,%s,%d,%s\n' \
    "$nucleus" "$decay_mode" "$sweep" "$ne" "$nu" "$nv" "$points" \
    "$wmin" "$ee" "$u" "$v" "$seconds" "$validation" "$exit_code" "$logfile" >> "$CSV"

  echo "    W_min=$wmin; validation=$validation; grid time=${seconds}s; exit=$exit_code"
  if [[ "$exit_code" -ne 0 && "$validation" != 'ERROR_NO_SCAN' ]]; then
    echo "    WARNING: Nonzero exit code. Inspect: $logfile" >&2
    tail -n 8 "$logfile" >&2 || true
  fi
  if [[ "$validation" == 'ERROR_NO_SCAN' ]]; then
    echo "    WARNING: No grid-scan report found. Inspect: $logfile" >&2
    tail -n 8 "$logfile" >&2 || true
  fi
}

for nucleus in 11C 3H; do
  if [[ "$nucleus" == '11C' ]]; then
    config="$PLUS_CONFIG"
    decay_mode='BetaPlusPolarised'
  else
    config="$MINUS_CONFIG"
    decay_mode='BetaMinusPolarised'
  fi

  if [[ "$MODE" == 'all' || "$MODE" == 'angular' ]]; then
    for n in "${ANGULAR_POINTS[@]}"; do
      run_one "$nucleus" "$decay_mode" "$config" angular "$ANGULAR_FIXED_NE" "$n" "$n"
    done
  fi
  if [[ "$MODE" == 'all' || "$MODE" == 'energy' ]]; then
    for n in "${ENERGY_POINTS[@]}"; do
      run_one "$nucleus" "$decay_mode" "$config" energy "$n" "$ENERGY_FIXED_NUV" "$ENERGY_FIXED_NUV"
    done
  fi
done

echo
echo "Finished. Results: $CSV"
echo "Logs: $OUT/logs"
echo "NOTE: scan_seconds is NA if the timing output was not added to ValidateTransitionOnce()."
