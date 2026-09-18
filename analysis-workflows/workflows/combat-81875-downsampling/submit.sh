#!/usr/bin/env bash
# Submit the 81875 downsampling and three-caller SV pipeline.
# Usage: ./submit.sh [--single] [--dry-run]

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../lib/load_config.sh"
export SVCFIT_CONFIG FORCE KEEP_INTERMEDIATE
cd "$SCRIPT_DIR"
mkdir -p logs

dry_run=false
single=false
for arg in "$@"; do
  case "$arg" in
    --dry-run) dry_run=true ;;
    --single) single=true ;;
    *) echo "ERROR: unknown option: $arg" >&2; exit 2 ;;
  esac
done

if $single; then array_spec=1; else array_spec="1-${DOWNSAMPLE_REPLICATES}"; fi
fake_id_file=$(mktemp)
echo 1000 > "$fake_id_file"
trap 'rm -f "$fake_id_file"' EXIT
submit() {
  if $dry_run; then
    local fake_id
    fake_id=$(( $(cat "$fake_id_file") + 1 ))
    echo "$fake_id" > "$fake_id_file"
    printf 'sbatch --parsable --array=%s' "$array_spec" >&2
    printf ' %q' "$@" >&2
    printf '\n' >&2
    echo "$fake_id"
  else
    sbatch --parsable --array="$array_spec" "$@"
  fi
}

down=$(submit 01_downsample_bam.sh)
delly=$(submit --dependency="afterok:${down}" 02_call_delly.sh)
gridss=$(submit --dependency="afterok:${down}" 02_call_gridss.sh)
manta=$(submit --dependency="afterok:${down}" 02_call_manta.sh)
g_delly=$(submit --dependency="afterok:${delly}" 03_genotype_delly.sh)
g_gridss=$(submit --dependency="afterok:${gridss}" 03_genotype_gridss.sh)
g_manta=$(submit --dependency="afterok:${manta}" 03_genotype_manta.sh)
survivor=$(submit --dependency="afterok:${g_delly}:${g_gridss}:${g_manta}" 04_merge_survivor.sh)
prepare=$(submit --dependency="afterok:${survivor}" 05_prepare_svcfit_vcf.sh)

printf 'downsample=%s delly=%s gridss=%s manta=%s survivor=%s prepare=%s\n' \
  "$down" "$delly" "$gridss" "$manta" "$survivor" "$prepare"
$dry_run && echo 'DRY RUN: no jobs submitted'
