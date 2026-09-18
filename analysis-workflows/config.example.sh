#!/usr/bin/env bash
# Central configuration for the preserved SVCFit analysis workflows.
#
# Copy this file to config.local.sh and edit values for the execution site.
# config.local.sh should remain untracked because it may contain protected paths.
# Workflow scripts should source lib/load_config.sh rather than embed site paths.

# Repository checkout and code provenance
export SVCFIT_REPO_ROOT="${SVCFIT_REPO_ROOT:-/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/SVCFit}"
# Fill this with the exact commit used for a run. Leave no production run unpinned.
export SVCFIT_EXPECTED_COMMIT="${SVCFIT_EXPECTED_COMMIT:-}"
export SVCFIT_R_SOURCE_DIR="${SVCFIT_R_SOURCE_DIR:-${SVCFIT_REPO_ROOT}/R}"

# Shared COMBAT project roots
export COMBAT_WGS_ROOT="${COMBAT_WGS_ROOT:-/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS}"
export SVCFIT_ANALYSIS_ROOT="${SVCFIT_ANALYSIS_ROOT:-${COMBAT_WGS_ROOT}/svcfit}"
export BAM_DIR="${BAM_DIR:-${COMBAT_WGS_ROOT}/bam_hg38}"
export FACETS_DIR="${FACETS_DIR:-${COMBAT_WGS_ROOT}/VY_COMBAT_WGS_000_analysis/facets}"
export GATK_GERMLINE_DIR="${GATK_GERMLINE_DIR:-${COMBAT_WGS_ROOT}/VY_COMBAT_WGS_000_analysis/HaplotypeCaller}"
export REFERENCE_FASTA="${REFERENCE_FASTA:-/projects/karchin-lab-hpc/hg38_reference/hg38.analysisSet.fa}"

# Cohort manifests
export PAIR_FILE="${PAIR_FILE:-${SVCFIT_ANALYSIS_ROOT}/pair.bed}"
export PURITY_FILE="${PURITY_FILE:-${SVCFIT_ANALYSIS_ROOT}/samp_pur.csv}"
export SAMPLE_MANIFEST="${SAMPLE_MANIFEST:-${SVCFIT_ANALYSIS_ROOT}/sv_call/svtyp/samples.txt}"
export EXCLUDED_PAIR="${EXCLUDED_PAIR:-84972__86626}"

# Input and output trees
export SV_CALL_ROOT="${SV_CALL_ROOT:-${SVCFIT_ANALYSIS_ROOT}/sv_call}"
export SURVIVOR_DIR="${SURVIVOR_DIR:-${SV_CALL_ROOT}/survivor}"
export ORIGINAL_SNP_DIR="${ORIGINAL_SNP_DIR:-${SVCFIT_ANALYSIS_ROOT}/snps}"
export ORIGINAL_OUTPUT_DIR="${ORIGINAL_OUTPUT_DIR:-${SVCFIT_ANALYSIS_ROOT}/output}"
export GERMLINE_RERUN_ROOT="${GERMLINE_RERUN_ROOT:-${SVCFIT_ANALYSIS_ROOT}/germline_rerun}"
export GERMLINE_HET_DIR="${GERMLINE_HET_DIR:-${GERMLINE_RERUN_ROOT}/germline}"
export GERMLINE_SNP_DIR="${GERMLINE_SNP_DIR:-${GERMLINE_RERUN_ROOT}/snps}"
export GERMLINE_OUTPUT_DIR="${GERMLINE_OUTPUT_DIR:-${GERMLINE_RERUN_ROOT}/output}"

# Conda and tool environments
export CONDA_SH="${CONDA_SH:-/home/yliu498/miniforge3/etc/profile.d/conda.sh}"
export ENV_ALIGN="${ENV_ALIGN:-/home/yliu498/miniforge3/envs/align}"
export ENV_VISOR="${ENV_VISOR:-/home/yliu498/miniforge3/envs/visor}"
export ENV_SVCFIT="${ENV_SVCFIT:-/home/yliu498/miniforge3/envs/SVCFit}"
export ENV_PYTHON="${ENV_PYTHON:-/home/yliu498/miniforge3/envs/py3}"
export ENV_DELLY="${ENV_DELLY:-/home/yliu498/miniforge3/envs/delly}"
export ENV_GRIDSS="${ENV_GRIDSS:-/home/yliu498/miniforge3/envs/gridss}"
export ENV_MANTA="${ENV_MANTA:-/home/yliu498/miniforge3/envs/manta}"
export ENV_SVTYPER="${ENV_SVTYPER:-/home/yliu498/miniforge3/envs/svtyp}"
export ENV_SURVIVOR="${ENV_SURVIVOR:-/home/yliu498/miniforge3/envs/survivor}"
export GRIDSS_JAR="${GRIDSS_JAR:-${ENV_GRIDSS}/share/gridss-2.13.2-6/gridss.jar}"

# SVCFit parameters shared by the COMBAT rerun
export SVCFIT_THRESHOLD="${SVCFIT_THRESHOLD:-0.1}"
export SVCFIT_FLANK_DEL="${SVCFIT_FLANK_DEL:-50}"
export SVCFIT_FLANK_SNP="${SVCFIT_FLANK_SNP:-500}"
export SVCFIT_FLANK_CNV="${SVCFIT_FLANK_CNV:-1000}"
export SVCFIT_QUAL_THRESHOLD="${SVCFIT_QUAL_THRESHOLD:-100}"
export SVCFIT_MIN_ALT="${SVCFIT_MIN_ALT:-2}"

# Sample-81875 downsampling analysis
export DOWNSAMPLE_TUMOR_SAMPLE="${DOWNSAMPLE_TUMOR_SAMPLE:-81875}"
export DOWNSAMPLE_NORMAL_SAMPLE="${DOWNSAMPLE_NORMAL_SAMPLE:-PBMC_10}"
export DOWNSAMPLE_ON_TREATMENT_SAMPLE="${DOWNSAMPLE_ON_TREATMENT_SAMPLE:-82780_recut}"
export DOWNSAMPLE_SOURCE_COVERAGE="${DOWNSAMPLE_SOURCE_COVERAGE:-36}"
export DOWNSAMPLE_TARGET_COVERAGE="${DOWNSAMPLE_TARGET_COVERAGE:-17}"
export DOWNSAMPLE_REPLICATES="${DOWNSAMPLE_REPLICATES:-20}"
export DOWNSAMPLE_SEED_BASE="${DOWNSAMPLE_SEED_BASE:-1000}"
export DOWNSAMPLE_OUTPUT_DIR="${DOWNSAMPLE_OUTPUT_DIR:-${BAM_DIR}/down_samp}"
export DOWNSAMPLE_TUMOR_PURITY="${DOWNSAMPLE_TUMOR_PURITY:-0.68}"
export DOWNSAMPLE_ON_TREATMENT_PURITY="${DOWNSAMPLE_ON_TREATMENT_PURITY:-0.37}"

# Runtime behavior
export FORCE="${FORCE:-0}"
export KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-0}"

