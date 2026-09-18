#!/usr/bin/env bash
#SBATCH --job-name=manta_svtyper
#SBATCH --output=logs/genotype_manta_%A_%a.out
#SBATCH --error=logs/genotype_manta_%A_%a.err
#SBATCH --time=20:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=3

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../lib/load_config.sh"
source "$CONDA_SH"
conda activate "$ENV_SVTYPER"

id=${BOOTSTRAP_ID:-${SLURM_ARRAY_TASK_ID:?Set BOOTSTRAP_ID or run as a SLURM array task}}
sample="boot${id}"
input="${SV_CALL_ROOT}/manta/${sample}/results/variants/p_${sample}.vcf"
bam="${DOWNSAMPLE_OUTPUT_DIR}/bootstrap_${id}/${DOWNSAMPLE_TUMOR_SAMPLE}_covdown_${DOWNSAMPLE_TARGET_COVERAGE}x_boot${id}.rg.bam"
out_dir="${SV_CALL_ROOT}/svtyp/${sample}"
tmp="${out_dir}/m${sample}.vcf"
output="${out_dir}/manta_${sample}.vcf"
for path in "$input" "$bam"; do [[ -f "$path" ]] || { echo "ERROR: missing input: $path" >&2; exit 1; }; done
mkdir -p "$out_dir"

awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {info=$8; gsub(/CIPOS=[^;]+;?/,"",info); gsub(/CIEND=[^;]+;?/,"",info); $8="CIPOS=-100,100;CIEND=-100,100;" info; print}' "$input" > "$tmp"
svtyper-sso --core "${SLURM_CPUS_PER_TASK:-3}" --batch_size 1000 --max_reads 1000 -i "$tmp" -B "$bam" > "$output"
rm -f "$tmp"
[[ -s "$output" ]] || { echo "ERROR: empty SVTyper output: $output" >&2; exit 1; }
