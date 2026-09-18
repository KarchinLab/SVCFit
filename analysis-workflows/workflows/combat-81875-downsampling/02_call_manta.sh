#!/usr/bin/env bash
#SBATCH --job-name=manta
#SBATCH --output=logs/manta_%A_%a.out
#SBATCH --error=logs/manta_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=15G
#SBATCH --time=20:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../lib/load_config.sh"
source "$CONDA_SH"
conda activate "$ENV_MANTA"

id=${BOOTSTRAP_ID:-${SLURM_ARRAY_TASK_ID:?Set BOOTSTRAP_ID or run as a SLURM array task}}
sample="boot${id}"
tumor_bam="${DOWNSAMPLE_OUTPUT_DIR}/bootstrap_${id}/${DOWNSAMPLE_TUMOR_SAMPLE}_covdown_${DOWNSAMPLE_TARGET_COVERAGE}x_boot${id}.rg.bam"
normal_bam="${BAM_DIR}/${DOWNSAMPLE_NORMAL_SAMPLE}_hg38_recal_sorted.bam"
out_dir="${SV_CALL_ROOT}/manta/${sample}"
for input in "$tumor_bam" "$normal_bam" "$REFERENCE_FASTA" "$MANTA_CONFIG" "$MANTA_CONVERT_INV"; do
  [[ -e "$input" ]] || { echo "ERROR: missing input/tool: $input" >&2; exit 1; }
done
[[ ! -e "$out_dir" || "$FORCE" == "1" ]] || { echo "ERROR: Manta run directory exists: $out_dir (set FORCE=1 to replace)" >&2; exit 1; }
[[ "$FORCE" != "1" ]] || rm -rf "$out_dir"
"$MANTA_CONFIG" --normalBam "$normal_bam" --tumorBam "$tumor_bam" \
  --referenceFasta "$REFERENCE_FASTA" --runDir "$out_dir"
python "${out_dir}/runWorkflow.py" -j "${SLURM_CPUS_PER_TASK:-3}"
python "$MANTA_CONVERT_INV" "$SAMTOOLS_BIN" "$REFERENCE_FASTA" \
  "${out_dir}/results/variants/somaticSV.vcf.gz" > "${out_dir}/results/variants/${sample}.vcf"
bcftools view -f PASS -Ov -o "${out_dir}/results/variants/p_${sample}.vcf" \
  "${out_dir}/results/variants/${sample}.vcf"
[[ -s "${out_dir}/results/variants/p_${sample}.vcf" ]] || { echo "ERROR: empty Manta output" >&2; exit 1; }
