#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

RUN_ID="${SERIALNULL_RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
export SERIALNULL_RUN_ID="${RUN_ID}"

if [[ -z "${R_BIN:-}" ]]; then
  if command -v Rscript >/dev/null 2>&1; then
    R_BIN="Rscript"
  elif [[ -x /usr/local/bin/Rscript ]]; then
    R_BIN="/usr/local/bin/Rscript"
  else
    R_BIN="Rscript"
  fi
fi
PYTHON_BIN="${PYTHON_BIN:-${REPO_ROOT}/analysis/.venv-graphcompass/bin/python}"

if ! command -v "${R_BIN}" >/dev/null 2>&1; then
  echo "[ERROR] Rscript not found (R_BIN=${R_BIN})" >&2
  exit 1
fi

if [[ ! -x "${PYTHON_BIN}" ]]; then
  echo "[ERROR] Python interpreter not found or not executable: ${PYTHON_BIN}" >&2
  echo "        Set PYTHON_BIN to a Python with graphcompass installed." >&2
  exit 1
fi

CACHE_BASE="${REPO_ROOT}/analysis/.cache-python"
mkdir -p "${CACHE_BASE}" "${CACHE_BASE}/matplotlib" "${CACHE_BASE}/fontconfig"
export XDG_CACHE_HOME="${CACHE_BASE}"
export MPLCONFIGDIR="${CACHE_BASE}/matplotlib"
export FONTCONFIG_PATH="${CACHE_BASE}/fontconfig"

cd "${REPO_ROOT}"

echo "[INFO] SERIALNULL_RUN_ID=${SERIALNULL_RUN_ID}"

echo "[INFO] Preflight: parse-check R scripts"
"${R_BIN}" -e "invisible(parse(file='analysis/serialnull/serialnull_utils.R')); invisible(parse(file='analysis/serialnull/01_prepare_manifest.R')); invisible(parse(file='analysis/serialnull/02_run_celledger.R')); invisible(parse(file='analysis/serialnull/04_summarize_timings.R')); cat('Preflight parse: OK\\n')"

echo "[INFO] Step 1/4: prepare shared sample manifest"
"${R_BIN}" analysis/serialnull/01_prepare_manifest.R

echo "[INFO] Step 2/4: run CellEdgeR analyses"
"${R_BIN}" analysis/serialnull/02_run_celledger.R

echo "[INFO] Step 3/4: run GraphCompass analyses"
"${PYTHON_BIN}" analysis/serialnull/03_run_graphcompass.py

echo "[INFO] Step 4/4: summarize timings"
"${R_BIN}" analysis/serialnull/04_summarize_timings.R

echo "[INFO] Done"
echo "       Manifest: analysis/serialnull/config/sample_splits.csv"
echo "       Timings log: analysis/serialnull/results/timings.csv"
echo "       Timing totals: analysis/serialnull/results/timings_totals_latest.csv"
echo "       CellEdgeR results: analysis/serialnull/results/celledger"
echo "       GraphCompass results: analysis/serialnull/results/graphcompass"
