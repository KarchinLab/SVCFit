# Review of `SVCFit-2024-2026/script`

Date: 18 September 2026

## Summary

The directory contains 21 shell scripts, 9 R scripts, and one useful session note. It is a mixture of two real analysis pipelines, older superseded scripts, and a few exploratory plotting files. All shell scripts pass `bash -n`, all R files parse, and the sample-81875 pipeline completes its built-in dry run. This checks syntax and orchestration only; no HPC jobs or analyses were run.

The strongest material to preserve is the germline-corrected COMBAT rerun and the documented sample-81875 downsampling pipeline. They should be kept as separate workflow directories with a shared configuration file, README, and pinned software versions. At present, 27 of 30 scripts contain hard-coded HPC or home-directory paths, and only 10 of 21 shell scripts enable strict shell error handling.

## Source availability

| Source | Status from this computer | Consequence |
|---|---|---|
| `SVCFit-2024-2026/script` | **Located.** All 32 files are readable locally and recorded in `legacy-script-directory.sha256`. | The legacy code can be reviewed and copied without changing its source directory. |
| `~/Documents/git-repos/SVCFit` | **Located and synchronized.** Local `main` and `origin/main` both point to commit `9c663a1c401b46243e7e8a84638b19725765e9ab` after fetching on 18 September 2026. | The configuration and preservation record can be committed to the current package repository. |
| `/projects/karchin-lab-hpc/COMBAT/COMBAT_WGS/...` | **Not available from this computer; protected-server location.** | BAM, VCF, FACETS, purity, pair, recovery, and completed output files referenced by the scripts could not be inspected or hashed. |
| `/projects/karchin-lab-hpc/hg38_reference/hg38.analysisSet.fa` | **Not available from this computer; protected-server location.** | The reference FASTA and its indexes could not be verified. |
| `/home/yliu498/miniforge3/...`, `~/manta/...`, and cluster executables | **Not available from this computer; server-specific software locations.** | Conda environments, GRIDSS, Manta, SVTyper, SURVIVOR, and their versions require verification on the cluster. |
| Mendeley dataset referenced by the package README | **Location documented; contents not downloaded for this review.** | Public deposited inputs and outputs were not compared with the legacy scripts. |
| European Genome-phenome Archive accession `EGAD00001001343` | **Protected source; not accessed.** | Controlled-access data were not inspected. |

These unavailable sources are configuration and validation dependencies, not missing code files. They should remain outside Git. The protected-server smoke tests in the migration plan are the point at which their existence, permissions, versions, and expected checksums should be recorded.

## Recommended disposition

| Disposition | Files | Reason |
|---|---|---|
| **Preserve as the current germline-correction workflow** | `combat_02b_filter_germline.sh`, `combat_03_snp_process.sh`, `submit_germline_rerun.sh`, `combat_svcfit_rerun.sh`, `run_svcfit_cluster_tree.R`, `get_sv_range.R`, `merge_sv2clust.R` | This is a coherent, documented normal-derived germline workflow. It protects the original outputs, has dependency chaining and idempotency checks, and directly supports the corrected COMBAT results. |
| **Preserve as an analysis-specific provenance workflow** | `run_pipeline.sh`, `run_down_samp_delly.sh`, `run_down_samp_gridss.sh`, `run_down_samp_manta.sh`, `del_svtyp.sh`, `gri_svtyp.sh`, `man_svtyp.sh`, `run_survivor.sh`, `proc_surv.sh`, `surv2svcfit.R`, `session_notes.md` | Together these record the 20-replicate, 17x sample-81875 downsampling and three-caller SV pipeline. The dependency graph dry-runs correctly. Keep the set together; individual scripts have little meaning outside it. |
| **Preserve after review or consolidation** | `run_svcfit_pair.sh`, `run_svcfit_cli.R`, `run_svcfit_boot.sh`, `extract_sv_breakpoint_depth.sh`, `old_down_cov.sh` | These contain useful general wrappers or analysis provenance, but overlap newer workflows or have inconsistencies described below. `old_down_cov.sh` is actually the safer and better-documented downsampler and should be renamed if retained. |
| **Archive as superseded, not runnable current code** | `run_snps.sh`, `svcfit.sh`, `get_read_for_surv.R`, `cluster.R` | `run_snps.sh` selects heterozygous sites from tumor calls and is superseded by the normal-derived germline correction. `svcfit.sh` is superseded and passes `True` to a flag-style option. `get_read_for_surv.R` duplicates the older deposited implementation and is superseded by `surv2svcfit.R`. `cluster.R` is a hard-coded exploratory predecessor of the package implementation. |
| **Discard or retain only in a clearly labeled scratch archive** | `AR_cnv.R`, `make_circus.R`, `run_circus.sh`, `.DS_Store` | `AR_cnv.R` is a one-off interactive plot with no saved output. The circos pair is not runnable as written: the shell wrapper supplies command-line options that the R script never parses, the R script lacks its library setup, refers to undefined `X` and `Y`, hard-codes sample 81875, and uses the superseded inversion/duplication color assignment. |

## Issues to fix before GitHub preservation

1. **Pin the executed SVCFit revision.** `run_svcfit_cluster_tree.R` and `run_svcfit_cli.R` source every R file from a fixed HPC checkout at runtime. The saved workflow should record a commit and load that checkout explicitly; otherwise the same script can produce different results after the package changes.

2. **Replace `down_cov_81875.sh` with the safer implementation.** The current 25-line file writes a BAM named `.rg.bam` but does not replace its read group. `old_down_cov.sh` does replace the read group, computes and validates the sampling fraction, records seeds, quotes paths, and optionally removes intermediates. The session notes describe this safer implementation, not the current short file.

3. **Correct `run_svcfit_boot.sh` before reuse.** Its comments say clustering and tree construction are skipped with `--stages svcfit`, but the command omits that option, so the default runs all three stages. The comments also say the same replicate is supplied twice, while the command supplies the bootstrap sample and `82780_recut`.

4. **Separate configuration from code.** Paths, conda environments, sample IDs, purities, array sizes, references, and tool locations should live in one checked configuration or manifest. The current scripts assume the exact `/projects/karchin-lab-hpc/...` layout and Yunzhou's environment paths.

5. **Add input and output validation to the older SV-calling scripts.** Much of the 81875 workflow lacks `set -euo pipefail`, quotes, uniqueness checks for `ls | grep`, and output checks. It remains valuable provenance, but should not be presented as portable production code until these are added.

6. **Document fixed-format assumptions.** `get_sv_range.R`, `surv2svcfit.R`, and several callers parse VCFs by fixed column counts and regular-expression positions. Those assumptions should be stated and checked at startup because a caller or VCF-version change can silently shift fields.

## Suggested GitHub layout

```text
analysis-workflows/
  README.md
  config.example.sh
  lib/
    load_config.sh
  manifests/
    combat_pairs.tsv
    sample_to_normal.tsv
    software_versions.tsv
  workflows/
    combat-germline-rerun/
      README.md
      submit.sh
      01_filter_normal_hets.sh
      02_pileup_tumor_at_normal_hets.sh
      03_run_svcfit_pairs.sh
      run_svcfit_cluster_tree.R
      get_sv_ranges.R
      merge_sv_to_cluster.R
    combat-81875-downsampling/
      README.md
      submit.sh
      01_downsample_bam.sh
      02_call_delly.sh
      02_call_gridss.sh
      02_call_manta.sh
      03_genotype_delly.sh
      03_genotype_gridss.sh
      03_genotype_manta.sh
      04_merge_survivor.sh
      05_prepare_svcfit_vcf.R
      06_run_svcfit.sh
    breakpoint-depth/
      README.md
      extract_sv_breakpoint_depth.sh
archive/
  superseded/
    README.md
  scratch/
    README.md
```

The top-level `config.example.sh` is the single site-configuration interface. It holds repository and data roots, manifests, conda environments, shared SVCFit parameters, and the 81875 downsampling settings. Each workflow sources `lib/load_config.sh`; the R entry points read the exported values with `Sys.getenv()`. A local `config.local.sh` is ignored by Git. Protected data and machine-specific paths therefore stay outside version control while the parameter names and expected inputs remain visible.

Keep cohort membership and sample-to-normal relationships in versioned TSV manifests rather than shell arrays or `ls | grep` searches. Keep software versions separate from paths so each completed run can record the SVCFit commit and external tool versions. Generated data, logs, BAMs, VCFs, RDS files, and credentials should remain outside the repository.

## Preservation and migration plan

The original directory remains at `SVCFit-2024-2026/script` and has not been modified. Its 32 files are recorded in `script_preservation/legacy-script-directory.sha256`, which provides an exact baseline for later comparisons.

The centralized configuration is in `config.example.sh`, with validation and optional SVCFit commit checking in `lib/load_config.sh`. These files do not alter the legacy scripts.

Migration should proceed in this order:

1. Create the repository structure above and commit the untouched checksum manifest.
2. Move copies of the germline-correction workflow into `workflows/combat-germline-rerun`, replace hard-coded paths with configuration variables, and verify output hashes or scientific summaries against the completed corrected run.
3. Move the 81875 workflow as a complete unit. Use the implementation currently named `old_down_cov.sh` as the basis for `01_downsample_bam.sh`; retain the short `down_cov_81875.sh` only in the superseded archive.
4. Resolve the `run_svcfit_boot.sh` stage mismatch before installing it as step 06. Record whether the intended product is per-replicate SVCFit output alone or a paired clustering/tree analysis.
5. Put `run_snps.sh`, `svcfit.sh`, `get_read_for_surv.R`, and `cluster.R` in `archive/superseded` with a README explaining their replacements. Put the AR and circos experiments in `archive/scratch`, or omit them after the legacy snapshot is safely retained.
6. Run syntax checks, the 81875 dry run, configuration validation, and small protected-server smoke tests. Only then label the migrated workflows as current.

No files in the source `script` directory were changed during this review or configuration work.
