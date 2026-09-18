# COMBAT germline-corrected SVCFit rerun

This workflow rebuilds the COMBAT SVCFit analysis using germline heterozygous SNP sites selected from each subject's PBMC normal. Tumor BAMs provide read counts at those fixed sites; they do not decide which sites exist.

## Steps

1. `01_filter_normal_hets.sh` filters each PBMC HaplotypeCaller VCF to biallelic heterozygous SNPs.
2. `02_pileup_tumor_at_normal_hets.sh` selects sites near SVs and measures full-BAM and discordant/supplementary tumor support.
3. `03_run_svcfit_pairs.sh` runs SVCFit for each longitudinal pair, followed by clustering and tree construction.
4. `submit.sh` sizes each SLURM array from the configured manifests and submits the three stages with `afterok` dependencies.

`get_sv_ranges.R`, `run_svcfit_cluster_tree.R`, and `merge_sv_to_cluster.R` are workflow helpers. The paired runner loads SVCFit R sources from the commit-checked `SVCFIT_R_SOURCE_DIR`.

The configured R environment must provide `optparse`, `dplyr`, `tidyr`, `stringr`, `purrr`, `GenomicRanges`, `readr`, `ggplot2`, `RColorBrewer`, and `reticulate`. The shell steps require `bcftools`, `bgzip`, `tabix`, and `samtools` on `PATH`; the clustering stage also requires the Python dependencies used by SVCFit.

## Required configuration

Copy `../../config.example.sh` to `../../config.local.sh`, edit the protected-server values, and set `SVCFIT_EXPECTED_COMMIT`. The sample manifest must contain `sample matched_normal` columns. The pair file must contain two sample columns. The purity file must contain sample and purity columns separated by tabs.

Submit with:

```bash
export SVCFIT_CONFIG=/absolute/path/to/analysis-workflows/config.local.sh
./submit.sh
```

Set `FORCE=1` to rebuild existing products. Before a full run, test one normal, sample, and pair interactively with `NORMAL_SAMPLE`, `SAMPLE`, and `PAIR` respectively.

The scripts have been checked locally for syntax and parsing. Protected inputs and scientific equivalence with the completed rerun must be verified on the cluster.
