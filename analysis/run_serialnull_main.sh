#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_SCRIPT="${SCRIPT_DIR}/serialnull/run_serialnull_pipeline.sh"

if [[ ! -x "${PIPELINE_SCRIPT}" ]]; then
  echo "[ERROR] Missing pipeline script: ${PIPELINE_SCRIPT}" >&2
  exit 1
fi

exec "${PIPELINE_SCRIPT}" "$@"
