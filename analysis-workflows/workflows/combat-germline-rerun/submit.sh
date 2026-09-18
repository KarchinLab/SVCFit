#!/usr/bin/env bash
# Submit normal filtering, tumor pileup, and paired SVCFit as an afterok chain.
# Usage: ./submit.sh [--dry-run]

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../lib/load_config.sh"
export SVCFIT_CONFIG FORCE KEEP_INTERMEDIATE
cd "$SCRIPT_DIR"
mkdir -p logs

dry_run=false
for arg in "$@"; do
  case "$arg" in
    --dry-run) dry_run=true ;;
    *) echo "ERROR: unknown option: $arg" >&2; exit 2 ;;
  esac
done

normal_count=$(awk 'NF >= 2 && $1 !~ /^#/ {print $2}' "$SAMPLE_MANIFEST" | sort -u | wc -l | tr -d ' ')
sample_count=$(awk -v excluded="$EXCLUDED_PAIR" 'NF >= 2 && $1 !~ /^#/ && $1"__"$2 != excluded {print $1; print $2}' "$PAIR_FILE" | awk '!seen[$0]++' | wc -l | tr -d ' ')
pair_count=$(awk -v excluded="$EXCLUDED_PAIR" 'NF >= 2 && $1 !~ /^#/ && $1"__"$2 != excluded {n++} END{print n+0}' "$PAIR_FILE")
((normal_count > 0 && sample_count > 0 && pair_count > 0)) || { echo "ERROR: one or more manifests are empty" >&2; exit 1; }

fake_id_file=$(mktemp)
echo 2000 > "$fake_id_file"
trap 'rm -f "$fake_id_file"' EXIT
submit() {
  if $dry_run; then
    local fake_id
    fake_id=$(( $(cat "$fake_id_file") + 1 ))
    echo "$fake_id" > "$fake_id_file"
    printf 'sbatch --parsable' >&2
    printf ' %q' "$@" >&2
    printf '\n' >&2
    echo "$fake_id"
  else
    sbatch --parsable "$@"
  fi
}

normal_job=$(submit --array="0-$((normal_count-1))" 01_filter_normal_hets.sh)
sample_job=$(submit --array="0-$((sample_count-1))" --dependency="afterok:${normal_job}" 02_pileup_tumor_at_normal_hets.sh)
pair_job=$(submit --array="0-$((pair_count-1))" --dependency="afterok:${sample_job}" 03_run_svcfit_pairs.sh)
printf 'normal_filter=%s tumor_pileup=%s svcfit_pairs=%s\n' "$normal_job" "$sample_job" "$pair_job"
$dry_run && echo 'DRY RUN: no jobs submitted'
