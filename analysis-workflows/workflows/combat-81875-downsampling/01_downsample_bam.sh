#!/usr/bin/env bash
#SBATCH --job-name=cov_down_81875
#SBATCH --output=logs/downsample_%A_%a.out
#SBATCH --error=logs/downsample_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../lib/load_config.sh"
source "$CONDA_SH"
conda activate "$ENV_ALIGN"

boot=${BOOTSTRAP_ID:-${SLURM_ARRAY_TASK_ID:?Set BOOTSTRAP_ID or run as a SLURM array task}}
threads=${SLURM_CPUS_PER_TASK:-4}
sample_name="${DOWNSAMPLE_TUMOR_SAMPLE}_covdown_${DOWNSAMPLE_TARGET_COVERAGE}x_boot${boot}"
tumor_bam="${BAM_DIR}/${DOWNSAMPLE_TUMOR_SAMPLE}_hg38_recal_sorted.bam"
out_dir="${DOWNSAMPLE_OUTPUT_DIR}/bootstrap_${boot}"
sub_bam="${out_dir}/${sample_name}.tumor_sub.bam"
final_bam="${out_dir}/${sample_name}.rg.bam"

[[ -f "$tumor_bam" ]] || { echo "ERROR: missing tumor BAM: $tumor_bam" >&2; exit 1; }
mkdir -p "$out_dir"
if [[ -s "$final_bam" && "$FORCE" != "1" ]]; then
  echo "[$sample_name] output exists; set FORCE=1 to rebuild: $final_bam"
  exit 0
fi
fraction=$(awk -v target="$DOWNSAMPLE_TARGET_COVERAGE" -v source="$DOWNSAMPLE_SOURCE_COVERAGE" \
  'BEGIN { if (source <= 0) exit 2; printf "%.6f", target/source }')
awk -v x="$fraction" 'BEGIN { if (x <= 0 || x > 1) { print "ERROR: downsampling fraction out of range: " x > "/dev/stderr"; exit 1 } }'
seed=$((DOWNSAMPLE_SEED_BASE + boot))
samtools_fraction="${seed}.${fraction#0.}"

echo "[$(date)] replicate=$boot sample=$sample_name fraction=$fraction seed=$seed"
samtools view -@ "$threads" -b -s "$samtools_fraction" "$tumor_bam" > "$sub_bam"
samtools index -@ "$threads" "$sub_bam"
samtools addreplacerg -@ "$threads" \
  -r "@RG\tID:${DOWNSAMPLE_TUMOR_SAMPLE}cov_boot${boot}\tSM:${sample_name}\tLB:lib1\tPL:ILLUMINA\tPU:unit1" \
  -m overwrite_all -o "$final_bam" "$sub_bam"
samtools index -@ "$threads" "$final_bam"
samtools quickcheck "$final_bam"
if [[ "$KEEP_INTERMEDIATE" == "0" ]]; then rm -f "$sub_bam" "${sub_bam}.bai"; fi
echo "[$(date)] wrote $final_bam"
