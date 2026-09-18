#!/usr/bin/env bash
#SBATCH --job-name=combat_svcfit
#SBATCH --output=logs/svcfit_pairs_%A_%a.out
#SBATCH --error=logs/svcfit_pairs_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../lib/load_config.sh"
source "$CONDA_SH"
conda activate "$ENV_VISOR"

if [[ -n "${PAIR:-}" ]]; then
  pre=${PAIR%%__*}; on=${PAIR#*__}
else
  index=${SLURM_ARRAY_TASK_ID:?Set PAIR or run as a SLURM array task}
  mapfile -t pairs < <(awk -v excluded="$EXCLUDED_PAIR" 'NF >= 2 && $1 !~ /^#/ && $1"__"$2 != excluded {print $1"\t"$2}' "$PAIR_FILE")
  IFS=$'\t' read -r pre on <<< "${pairs[$index]:-}"
fi
[[ -n "${pre:-}" && -n "${on:-}" ]] || { echo "ERROR: no pair selected" >&2; exit 1; }
pair="${pre}__${on}"
out_dir="${GERMLINE_OUTPUT_DIR}/${pair}"
mkdir -p "$out_dir"
if [[ -s "${out_dir}/clustering/cluster_centroids.csv" && "$FORCE" != "1" ]]; then
  echo "[$pair] clustering exists; set FORCE=1 to rebuild"
  exit 0
fi

lookup_normal() { awk -v s="$1" 'NF >= 2 && $1 == s {print $2; exit}' "$SAMPLE_MANIFEST"; }
lookup_purity() { awk -v s="$1" 'BEGIN{FS="\t"} $1 == s {print $2; exit}' "$PURITY_FILE"; }
resolve_inputs() {
  local sample=$1 normal
  normal=$(lookup_normal "$sample")
  [[ -n "$normal" ]] || { echo "ERROR: no normal for $sample" >&2; return 1; }
  R_HET="${GERMLINE_SNP_DIR}/het_near_sv_${sample}.vcf"
  R_ON="${GERMLINE_SNP_DIR}/het_on_sv_${sample}.vcf"
  R_SV="${SURVIVOR_DIR}/${sample}/svcfit_${sample}.vcf"
  R_CNV="${FACETS_DIR}/${sample}_vs_${normal}_FACETS_cncf.txt"
  R_PURITY=$(lookup_purity "$sample")
  for path in "$R_HET" "$R_ON" "$R_SV" "$R_CNV"; do [[ -f "$path" ]] || { echo "ERROR [$sample]: missing $path" >&2; return 1; }; done
  [[ -n "$R_PURITY" && "$R_PURITY" != "NA" ]] || { echo "ERROR [$sample]: missing purity" >&2; return 1; }
}

resolve_inputs "$pre"; het1=$R_HET; on1=$R_ON; sv1=$R_SV; cnv1=$R_CNV; purity1=$R_PURITY
resolve_inputs "$on"; het2=$R_HET; on2=$R_ON; sv2=$R_SV; cnv2=$R_CNV; purity2=$R_PURITY
Rscript "${SCRIPT_DIR}/run_svcfit_cluster_tree.R" \
  --het_t1 "$het1" --on_t1 "$on1" --sv_t1 "$sv1" --cnv_t1 "$cnv1" --sample_t1 "$pre" --purity_t1 "$purity1" \
  --het_t2 "$het2" --on_t2 "$on2" --sv_t2 "$sv2" --cnv_t2 "$cnv2" --sample_t2 "$on" --purity_t2 "$purity2" \
  --thresh "$SVCFIT_THRESHOLD" --out_dir "$out_dir" --exper "$pair" --tum_only --python_env "$ENV_PYTHON"
echo "[$pair] complete"
