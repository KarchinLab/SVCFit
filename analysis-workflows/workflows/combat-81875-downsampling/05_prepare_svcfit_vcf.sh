#!/usr/bin/env bash
#SBATCH --job-name=prepare_svcfit_vcf
#SBATCH --output=logs/prepare_svcfit_%A_%a.out
#SBATCH --error=logs/prepare_svcfit_%A_%a.err
#SBATCH --time=3:00:00
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../lib/load_config.sh"
source "$CONDA_SH"
conda activate "$ENV_SVCFIT"

id=${BOOTSTRAP_ID:-${SLURM_ARRAY_TASK_ID:?Set BOOTSTRAP_ID or run as a SLURM array task}}
sample="boot${id}"
svtyp_dir="${SV_CALL_ROOT}/svtyp/${sample}"
out_dir="${SURVIVOR_DIR}/${sample}"
survivor_vcf="${out_dir}/${sample}_wbnd.vcf"
delly="${svtyp_dir}/delly_${sample}.vcf"
gridss="${svtyp_dir}/gridss_${sample}.vcf"
manta="${svtyp_dir}/manta_${sample}.vcf"
recovery="${SURVIVOR_DIR}/recover_sv/${sample}"
output="${out_dir}/svcfit_${sample}.vcf"
for path in "$survivor_vcf" "$delly" "$gridss" "$manta"; do [[ -s "$path" ]] || { echo "ERROR: missing or empty input: $path" >&2; exit 1; }; done

args=(-s "$survivor_vcf" -d "$delly" -g "$gridss" -m "$manta" -o "$output")
if [[ -n "$DOWNSAMPLE_REFERENCE_DELLY_VCF" && -s "$recovery" ]]; then
  [[ -s "$DOWNSAMPLE_REFERENCE_DELLY_VCF" ]] || { echo "ERROR: missing recovery reference: $DOWNSAMPLE_REFERENCE_DELLY_VCF" >&2; exit 1; }
  args+=(-r "$DOWNSAMPLE_REFERENCE_DELLY_VCF" -j "$recovery")
fi
Rscript "${SCRIPT_DIR}/05_prepare_svcfit_vcf.R" "${args[@]}"
[[ -s "$output" ]] || { echo "ERROR: empty prepared SVCFit VCF: $output" >&2; exit 1; }
