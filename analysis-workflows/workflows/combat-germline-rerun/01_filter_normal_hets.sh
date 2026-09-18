#!/usr/bin/env bash
#SBATCH --job-name=combat_germfilter
#SBATCH --output=logs/filter_normal_hets_%A_%a.out
#SBATCH --error=logs/filter_normal_hets_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../lib/load_config.sh"
source "$CONDA_SH"
conda activate "$ENV_VISOR"

if [[ -n "${NORMAL_SAMPLE:-}" ]]; then
  normal=$NORMAL_SAMPLE
else
  index=${SLURM_ARRAY_TASK_ID:?Set NORMAL_SAMPLE or run as a SLURM array task}
  mapfile -t normals < <(awk 'NF >= 2 && $1 !~ /^#/ {print $2}' "$SAMPLE_MANIFEST" | sort -u)
  normal=${normals[$index]:-}
fi
[[ -n "${normal:-}" ]] || { echo "ERROR: no normal sample selected" >&2; exit 1; }
[[ "$normal" == PBMC_* ]] || { echo "ERROR: expected PBMC normal, got: $normal" >&2; exit 1; }

provided="${GATK_GERMLINE_DIR}/${normal}_haplotypecaller.vcf"
output="${GERMLINE_HET_DIR}/${normal}_germline_het.vcf.gz"
[[ -f "$provided" ]] || { echo "ERROR: missing normal VCF: $provided" >&2; exit 1; }
mkdir -p "$GERMLINE_HET_DIR"
if [[ -s "$output" && "$FORCE" != "1" ]]; then
  echo "[$normal] output exists; set FORCE=1 to rebuild: $output"
  exit 0
fi

bcftools view -v snps -m2 -M2 -g het -Oz -o "$output" "$provided"
tabix -f -p vcf "$output"
echo "[$normal] $(bcftools view -H "$output" | wc -l | tr -d ' ') biallelic heterozygous sites"
