#!/usr/bin/env bash
#SBATCH --job-name=combat_snpproc
#SBATCH --output=logs/pileup_tumor_%A_%a.out
#SBATCH --error=logs/pileup_tumor_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=4:00:00

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/../../lib/load_config.sh"
source "$CONDA_SH"
conda activate "$ENV_VISOR"

if [[ -n "${SAMPLE:-}" ]]; then
  sample=$SAMPLE
else
  index=${SLURM_ARRAY_TASK_ID:?Set SAMPLE or run as a SLURM array task}
  mapfile -t samples < <(awk -v excluded="$EXCLUDED_PAIR" 'NF >= 2 && $1 !~ /^#/ && $1"__"$2 != excluded {print $1; print $2}' "$PAIR_FILE" | awk '!seen[$0]++')
  sample=${samples[$index]:-}
fi
[[ -n "${sample:-}" ]] || { echo "ERROR: no tumor sample selected" >&2; exit 1; }
normal=$(awk -v s="$sample" 'NF >= 2 && $1 == s {print $2; exit}' "$SAMPLE_MANIFEST")
[[ -n "$normal" ]] || { echo "ERROR: no matched normal for $sample in $SAMPLE_MANIFEST" >&2; exit 1; }

tumor_bam="${BAM_DIR}/${sample}_hg38_recal_sorted.bam"
germline_het="${GERMLINE_HET_DIR}/${normal}_germline_het.vcf.gz"
sv_vcf="${SURVIVOR_DIR}/${sample}/svcfit_${sample}.vcf"
out_near="${GERMLINE_SNP_DIR}/het_near_sv_${sample}.vcf"
out_on="${GERMLINE_SNP_DIR}/het_on_sv_${sample}.vcf"
for path in "$tumor_bam" "$germline_het" "$sv_vcf" "$REFERENCE_FASTA"; do
  [[ -f "$path" ]] || { echo "ERROR [$sample]: missing input: $path" >&2; exit 1; }
done
mkdir -p "$GERMLINE_SNP_DIR"
if [[ -s "$out_on" && "$FORCE" != "1" ]]; then
  echo "[$sample] output exists; set FORCE=1 to rebuild"
  exit 0
fi

Rscript "${SCRIPT_DIR}/get_sv_ranges.R" -s "$sample" -o "$GERMLINE_SNP_DIR" -v "$sv_vcf"
regions="${GERMLINE_SNP_DIR}/${sample}.bed"
germ_near="${GERMLINE_SNP_DIR}/germ_near_${sample}.vcf.gz"
alleles="${GERMLINE_SNP_DIR}/germ_alleles_${sample}.tsv"
sites="${GERMLINE_SNP_DIR}/germ_sites_${sample}.bed"

bcftools view -R "$regions" "$germline_het" -Oz -o "$germ_near"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$germ_near" > "$alleles"
[[ -s "$alleles" ]] || { echo "ERROR [$sample]: no germline sites overlap SV regions" >&2; exit 1; }
awk 'BEGIN{OFS="\t"}{print $1,$2,$3","$4}' "$alleles" | bgzip -c > "${alleles}.gz"
tabix -f -s1 -b2 -e2 "${alleles}.gz"
awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' "$alleles" > "$sites"

raw_near="${GERMLINE_SNP_DIR}/tumor_at_germline_${sample}.vcf"
bcftools mpileup -f "$REFERENCE_FASTA" -a FORMAT/AD,FORMAT/DP -A -q15 -Q20 -R "$sites" -Ou "$tumor_bam" \
  | bcftools call -m -C alleles -T "${alleles}.gz" -Ov -o "$raw_near"
{
  grep '^##' "$raw_near"
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' "$sample"
  awk 'BEGIN{OFS="\t"} NR==FNR{k=$1":"$2;ref[k]=$3;alt[k]=$4;next} /^#/{next} {k=$1":"$2;if(!(k in ref))next;split($10,s,":");dp=s[3]+0;n=split(s[2],a,",");r=a[1]+0;v=(n>=2?a[2]+0:0);print $1,$2,".",ref[k],alt[k],".",".","DP="dp,"GT:AD:DP","0/1:"r","v":"dp}' "$alleles" "$raw_near"
} > "$out_near"
if grep -v '^#' "$out_near" | awk -F'\t' '$10 !~ /^[^:]+:[0-9]+,[0-9]+:[0-9]+$/' | grep -q .; then
  echo "ERROR [$sample]: output contains fields incompatible with parse_het_snps" >&2; exit 1
fi

supplementary="${GERMLINE_SNP_DIR}/sup_${sample}.bam"
samtools view -f 1 -F 2 -b "$tumor_bam" > "$supplementary"
samtools index "$supplementary"
bcftools mpileup -f "$REFERENCE_FASTA" -a DP -A -R "$sites" "$supplementary" -Ov > "$out_on"
rm -f "$supplementary" "${supplementary}.bai"
echo "[$sample] wrote $out_near and $out_on"
