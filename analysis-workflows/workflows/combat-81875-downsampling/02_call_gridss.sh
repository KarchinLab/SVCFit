#!/usr/bin/env bash
#SBATCH --job-name=gridss
#SBATCH --output=logs/gridss_%A_%a.out
#SBATCH --error=logs/gridss_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=45G
#SBATCH --time=20:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../lib/load_config.sh"
source "$CONDA_SH"
conda activate "$ENV_GRIDSS"

id=${BOOTSTRAP_ID:-${SLURM_ARRAY_TASK_ID:?Set BOOTSTRAP_ID or run as a SLURM array task}}
sample="boot${id}"
tumor_bam="${DOWNSAMPLE_OUTPUT_DIR}/bootstrap_${id}/${DOWNSAMPLE_TUMOR_SAMPLE}_covdown_${DOWNSAMPLE_TARGET_COVERAGE}x_boot${id}.rg.bam"
normal_bam="${BAM_DIR}/${DOWNSAMPLE_NORMAL_SAMPLE}_hg38_recal_sorted.bam"
out_dir="${SV_CALL_ROOT}/gridss/${sample}"
for input in "$tumor_bam" "$normal_bam" "$REFERENCE_FASTA" "$GRIDSS_JAR"; do
  [[ -f "$input" ]] || { echo "ERROR: missing input: $input" >&2; exit 1; }
done
mkdir -p "$out_dir"
gridss -r "$REFERENCE_FASTA" -j "$GRIDSS_JAR" -o "${out_dir}/${sample}.vcf" \
  -t "${SLURM_CPUS_PER_TASK:-3}" -w "$out_dir" "$normal_bam" "$tumor_bam"
bcftools view -f PASS -Ov -o "${out_dir}/p_${sample}.vcf" "${out_dir}/${sample}.vcf"
[[ -s "${out_dir}/p_${sample}.vcf" ]] || { echo "ERROR: empty GRIDSS output" >&2; exit 1; }
