#!/usr/bin/env Rscript
###############################################################################
# run_svcfit_cluster_tree.R
#
# Paired-sample pipeline: SVCFit inference for two timepoints (t1 and t2),
# DP-GMM clustering, and tree building. Adapted from the original single-sample
# SVCFit wrapper and extended to mirror the longitudinal pipeline:
#
#   Stage 1 (svcfit)        : run_svcfit() for t1 and t2 independently.
#                             BEDs written to <out_dir>/COMBAT/SVCFit_output/
#                             so cluster_data() can locate them via data_dir.
#   Stage 2+3 (cluster+tree): build_trees(run_clustering=TRUE, run_tree=TRUE)
#
# Example:
#   Rscript run_svcfit_cluster_tree.R \
#     --het_t1 t1.het.vcf --on_t1 t1.onsv.vcf --sv_t1 t1.sv.vcf --cnv_t1 t1.bed \
#     --het_t2 t2.het.vcf --on_t2 t2.onsv.vcf --sv_t2 t2.sv.vcf --cnv_t2 t2.bed \
#     --sample_t1 caseA_t1 --sample_t2 caseA_t2 \
#     --purity_t1 0.6     --purity_t2 0.6 \
#     --out_dir /path/to/out --exper exp1
###############################################################################
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(GenomicRanges)
  library(readr)
  library(ggplot2)
  library(RColorBrewer)
  library(reticulate)
})

option_list <- list(
  # ---- Per-timepoint inputs (t1) ----
  make_option(c("--het_t1"),    type = "character", help = "t1 heterozygous SNP VCF"),
  make_option(c("--on_t1"),     type = "character", help = "t1 SNPs-on-SV VCF"),
  make_option(c("--sv_t1"),     type = "character", help = "t1 SV VCF"),
  make_option(c("--cnv_t1"),    type = "character", help = "t1 CNV file (FACETS tab-delimited)"),
  make_option(c("--sample_t1"), type = "character", help = "t1 sample name"),
  make_option(c("--purity_t1"), type = "double",    help = "t1 tumor purity (0-1)"),

  # ---- Per-timepoint inputs (t2) ----
  make_option(c("--het_t2"),    type = "character", help = "t2 heterozygous SNP VCF"),
  make_option(c("--on_t2"),     type = "character", help = "t2 SNPs-on-SV VCF"),
  make_option(c("--sv_t2"),     type = "character", help = "t2 SV VCF"),
  make_option(c("--cnv_t2"),    type = "character", help = "t2 CNV file (FACETS tab-delimited)"),
  make_option(c("--sample_t2"), type = "character", help = "t2 sample name"),
  make_option(c("--purity_t2"), type = "double",    help = "t2 tumor purity (0-1)"),

  # ---- Common SVCFit options (inherited from original script) ----
  make_option(c("-t", "--thresh"),  type = "double",   default = 0.1,
              help = "SVCFit CNV-order decision threshold [default: %default]"),
  make_option(c("-o", "--out_dir"), type = "character",
              help = "Base output directory (BEDs, clustering, tree go here)"),
  make_option(c("-e", "--exper"),   type = "character", default = "exp1",
              help = "Experiment identifier [default: %default]"),
  make_option(c("--tum_only"),      action = "store_true", default = FALSE,
              help = "SVs were called from tumor-only BAM [default: %default]"),
  make_option(c("--flank_del"),     type = "integer",  default = 50L,
              help = "Max distance (bp) for DEL/BND overlap [default: %default]"),
  make_option(c("--flank_snp"),     type = "integer",  default = 500L,
              help = "Flanking window (bp) for mapping SNPs to SVs [default: %default]"),
  make_option(c("--flank_cnv"),     type = "integer",  default = 1000L,
              help = "Flanking window (bp) for background CNV assignment [default: %default]"),
  make_option(c("--QUAL_thresh"),   type = "integer",  default = 100L,
              help = "Minimum SV quality score [default: %default]"),
  make_option(c("--min_alt"),       type = "integer",  default = 2L,
              help = "Minimum SV-supporting read count [default: %default]"),

  # ---- Stage control ----
  make_option(c("--stages"), type = "character", default = "svcfit,cluster,tree",
              help = "Comma-separated stages to run: svcfit,cluster,tree [default: %default]"),

  # ---- Clustering / tree options (from longitudinal pipeline) ----
  make_option(c("--ccf_floor"),                 type = "numeric", default = 0.1,
              help = "Zero CCF values below this threshold before clustering [default: %default]"),
  make_option(c("--concentration"),             type = "numeric", default = 1,
              help = "DP-GMM concentration parameter [default: %default]"),
  make_option(c("--min_dist"),                  type = "numeric", default = 0.2,
              help = "Minimum cluster distance for merging [default: %default]"),
  make_option(c("--lineage_precedence_thresh"), type = "numeric", default = 0.2,
              help = "Lineage precedence threshold for tree [default: %default]"),
  make_option(c("--sum_filter_thresh"),         type = "numeric", default = 0.2,
              help = "Sum-rule filter threshold for tree [default: %default]"),
  make_option(c("--linear_penalty"),            type = "numeric", default = 0.3,
              help = "Linear-tree penalty [default: %default]"),
  make_option(c("--no_dedup"),                  action = "store_true", default = FALSE,
              help = "Disable SV deduplication across samples [default: %default]"),

  # ---- Python env for DP-GMM clustering ----
  make_option(c("--python_env"), type = "character",
              default = Sys.getenv("ENV_PYTHON", unset = ""),
              help = "Conda env with Python deps for DP-GMM clustering [default: %default]")
)

parser <- OptionParser(option_list = option_list, add_help_option = TRUE)
opts   <- parse_args(parser)

required <- c(
  "het_t1", "on_t1", "sv_t1", "cnv_t1", "sample_t1", "purity_t1",
  "het_t2", "on_t2", "sv_t2", "cnv_t2", "sample_t2", "purity_t2", "out_dir"
)
missing <- required[vapply(required, function(x) is.null(opts[[x]]) ||
  (is.character(opts[[x]]) && !nzchar(opts[[x]])), logical(1))]
if (length(missing)) stop("Missing required options: ", paste(missing, collapse = ", "))
input_names <- c("het_t1", "on_t1", "sv_t1", "cnv_t1", "het_t2", "on_t2", "sv_t2", "cnv_t2")
missing_files <- input_names[!vapply(input_names, function(x) file.exists(opts[[x]]), logical(1))]
if (length(missing_files)) stop("Input files do not exist: ", paste(missing_files, collapse = ", "))

stages   <- trimws(strsplit(opts$stages, ",")[[1]])
unknown_stages <- setdiff(stages, c("svcfit", "cluster", "tree"))
if (length(unknown_stages)) stop("Unknown stages: ", paste(unknown_stages, collapse = ", "))
run_tree <- "tree" %in% stages

# Initialize Python only if clustering/tree stages will run
if (any(c("cluster", "tree") %in% stages)) {
  if (!nzchar(opts$python_env)) stop("--python_env or ENV_PYTHON is required for clustering/tree stages")
  use_condaenv(opts$python_env, required = TRUE)
}

# Source the pinned SVCFit checkout selected by the central configuration.
folder_path <- Sys.getenv("SVCFIT_R_SOURCE_DIR", unset = "")
if (!nzchar(folder_path) || !dir.exists(folder_path)) {
  stop("SVCFIT_R_SOURCE_DIR must name a readable SVCFit R source directory")
}
for (f in list.files(path = folder_path, pattern = "\\.[Rr]$", full.names = TRUE)) {
  source(f)
}

chr_list <- c(paste0("chr", 1:22), "chrX", "chrY")

# ---- Output directory layout ----
# cluster_data() expects BEDs at <data_dir>/COMBAT/SVCFit_output/<sample>.bed
if (is.null(opts$out_dir)) stop("--out_dir is required")
svcfit_dir <- opts$out_dir
combat_dir <- file.path(svcfit_dir, "COMBAT", "SVCFit_output")
dir.create(combat_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# Stage 1: SVCFit inference for t1 and t2
###############################################################################
run_one_svcfit <- function(tp_label, p_het, p_onsv, p_sv, p_cnv, samp_name) {
  message(sprintf("\n=== SVCFit: %s (%s) ===", samp_name, tp_label))
  result <- run_svcfit(
    p_het      = p_het,
    p_onsv     = p_onsv,
    p_sv       = p_sv,
    p_cnv      = p_cnv,
    samp       = samp_name,
    exper      = opts$exper,
    chr_lst    = chr_list,
    thresh     = opts$thresh,
    tum_only   = opts$tum_only,
    flank_del  = opts$flank_del,
    flank_snp  = opts$flank_snp,
    flank_cnv  = opts$flank_cnv,
    QUAL_thresh = opts$QUAL_thresh,
    min_alt    = opts$min_alt
  )
  # Write to both COMBAT layout (required by cluster_data) and flat dir
  bed_combat <- file.path(combat_dir, sprintf("%s.bed", samp_name))
  bed_flat   <- file.path(svcfit_dir, sprintf("%s.bed", samp_name))
  write_delim(result$svcf, bed_combat, delim = "\t", quote = "none", col_names = TRUE)
  write_delim(result$svcf, bed_flat,   delim = "\t", quote = "none", col_names = TRUE)
  message(sprintf("Written %d SVs to %s", nrow(result$svcf), bed_combat))
  invisible(result)
}

if ("svcfit" %in% stages) {
  run_one_svcfit("t1", opts$het_t1, opts$on_t1, opts$sv_t1, opts$cnv_t1, opts$sample_t1)
  run_one_svcfit("t2", opts$het_t2, opts$on_t2, opts$sv_t2, opts$cnv_t2, opts$sample_t2)
}

###############################################################################
# Stage 2 + 3: Clustering and tree via build_trees()
###############################################################################
if (any(c("cluster", "tree") %in% stages)) {

  # pair_path: no-header TSV with pre_BAT and on_BAT columns
  pair_path <- file.path(svcfit_dir, "pair_path.txt")
  write.table(
    data.frame(pre_BAT = opts$sample_t1, on_BAT = opts$sample_t2),
    pair_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
  )

  # pur_path: TSV with 'sample' and 'purity' columns
  pur_path <- file.path(svcfit_dir, "pur_path.txt")
  write.table(
    data.frame(sample = c(opts$sample_t1, opts$sample_t2),
               purity = c(opts$purity_t1, opts$purity_t2)),
    pur_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )

  message(sprintf("\n=== Clustering%s: %s + %s ===",
                  if (run_tree) " + tree" else "",
                  opts$sample_t1, opts$sample_t2))

  full_result <- build_trees(
    run_clustering            = TRUE,
    pair_path                 = pair_path,
    pur_path                  = pur_path,
    data_dir                  = svcfit_dir,
    pair_num                  = 1L,
    deduplicate               = !opts$no_dedup,
    concentration             = opts$concentration,
    min_dist                  = opts$min_dist,
    run_tree                  = run_tree,
    lineage_precedence_thresh = opts$lineage_precedence_thresh,
    sum_filter_thresh         = opts$sum_filter_thresh,
    ccf_floor                 = opts$ccf_floor,
    linear_penalty            = opts$linear_penalty
  )

  # ---- Save clustering outputs ----
  cluster_dir <- file.path(svcfit_dir, "clustering")
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
  cluster_result <- full_result[[1]][[1]]   # full cluster-assignment table
  sv2cluster     <- full_result[[1]][[2]]   # SV-to-cluster mapping
  clones         <- full_result[[1]][[3]]   # clone CCF table for the pair
  write_csv(cluster_result, file.path(cluster_dir, "cluster_result.csv"))
  write_csv(sv2cluster,     file.path(cluster_dir, "sv2cluster.csv"))
  write_csv(clones,         file.path(cluster_dir, "cluster_centroids.csv"))
  message(sprintf("Clustering output saved to %s", cluster_dir))

  # ---- Save tree outputs ----
  if (run_tree && !is.null(full_result[[2]])) {
    tree_dir <- file.path(svcfit_dir, "tree")
    dir.create(tree_dir, recursive = TRUE, showWarnings = FALSE)
    best_tree <- full_result[[2]][[1]]
    topo_str  <- paste(apply(best_tree, 1,
                             function(r) paste(r["parent"], r["child"], sep = "->")),
                       collapse = "; ")
    write_csv(best_tree, file.path(tree_dir, "tree_edges.csv"))
    writeLines(topo_str, file.path(tree_dir, "topology.txt"))
    saveRDS(full_result[[2]], file.path(tree_dir, "tree_result.rds"))
    message(sprintf("Tree output saved to %s", tree_dir))
    message(sprintf("Topology: %s", topo_str))
  }
}

message(sprintf("\n=== Pipeline complete: %s + %s ===",
                opts$sample_t1, opts$sample_t2))
