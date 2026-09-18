#!/usr/bin/env bash
#SBATCH --job-name=survivor
#SBATCH --output=logs/survivor_%A_%a.out
#SBATCH --error=logs/survivor_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=10:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../lib/load_config.sh"
source "$CONDA_SH"
conda activate "$ENV_SURVIVOR"

id=${BOOTSTRAP_ID:-${SLURM_ARRAY_TASK_ID:?Set BOOTSTRAP_ID or run as a SLURM array task}}
sample="boot${id}"
svtyp_dir="${SV_CALL_ROOT}/svtyp/${sample}"
out_dir="${SURVIVOR_DIR}/${sample}"
inputs=("${svtyp_dir}/delly_${sample}.vcf" "${svtyp_dir}/gridss_${sample}.vcf" "${svtyp_dir}/manta_${sample}.vcf")
for path in "${inputs[@]}"; do [[ -s "$path" ]] || { echo "ERROR: missing or empty caller VCF: $path" >&2; exit 1; }; done
mkdir -p "$out_dir"
list="${out_dir}/vcf_${sample}.txt"
printf '%s\n' "${inputs[@]}" > "$list"
SURVIVOR merge "$list" 300 2 1 1 0 50 "${out_dir}/${sample}.vcf"
sed 's/SVTYPE=TRA/SVTYPE=BND/g' "${out_dir}/${sample}.vcf" > "${out_dir}/${sample}_wbnd.vcf"
[[ -s "${out_dir}/${sample}_wbnd.vcf" ]] || { echo "ERROR: empty SURVIVOR output" >&2; exit 1; }
