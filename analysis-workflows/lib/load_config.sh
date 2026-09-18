#!/usr/bin/env bash
# Source from workflow scripts:
#   SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
#   source "${SCRIPT_DIR}/../../lib/load_config.sh"

set -euo pipefail

PRESERVATION_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SVCFIT_CONFIG=${SVCFIT_CONFIG:-${PRESERVATION_ROOT}/config.local.sh}

if [[ ! -f "$SVCFIT_CONFIG" ]]; then
    echo "ERROR: configuration file not found: $SVCFIT_CONFIG" >&2
    echo "Copy ${PRESERVATION_ROOT}/config.example.sh to config.local.sh and edit it." >&2
    return 1 2>/dev/null || exit 1
fi

# shellcheck source=/dev/null
source "$SVCFIT_CONFIG"

required_vars=(
    SVCFIT_REPO_ROOT SVCFIT_R_SOURCE_DIR COMBAT_WGS_ROOT SVCFIT_ANALYSIS_ROOT
    BAM_DIR FACETS_DIR GATK_GERMLINE_DIR REFERENCE_FASTA PAIR_FILE PURITY_FILE
    SV_CALL_ROOT SURVIVOR_DIR CONDA_SH ENV_VISOR ENV_SVCFIT ENV_PYTHON
)

for name in "${required_vars[@]}"; do
    if [[ -z "${!name:-}" ]]; then
        echo "ERROR: required configuration variable is empty: $name" >&2
        return 1 2>/dev/null || exit 1
    fi
done

if [[ -n "${SVCFIT_EXPECTED_COMMIT:-}" ]] && command -v git >/dev/null 2>&1; then
    actual_commit=$(git -C "$SVCFIT_REPO_ROOT" rev-parse HEAD 2>/dev/null || true)
    if [[ "$actual_commit" != "$SVCFIT_EXPECTED_COMMIT" ]]; then
        echo "ERROR: SVCFit checkout is ${actual_commit:-unreadable}; expected $SVCFIT_EXPECTED_COMMIT" >&2
        return 1 2>/dev/null || exit 1
    fi
fi

