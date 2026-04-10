#' Run the complete SVCFit pipeline
#'
#' Master function that runs SVCF inference and, optionally, DP-GMM clustering
#' and phylogenetic tree building.
#'
#' @section Pipeline stages:
#' \enumerate{
#'   \item \strong{SVCF inference} (always runs): \code{extract_info} →
#'         \code{characterize_sv} → \code{calc_svcf}.
#'   \item \strong{Clustering} (when \code{run_clustering = TRUE}): runs
#'         \code{cluster_data} on paired pre/on-treatment samples using a
#'         Dirichlet-process Gaussian mixture model (requires Python sklearn
#'         via reticulate).
#'   \item \strong{Tree building} (when \code{run_tree = TRUE}): infers the
#'         best-scoring phylogenetic spanning tree from cluster CCFs via
#'         \code{build_tree}. Requires \code{run_clustering = TRUE}.
#' }
#'
#' @param p_het Character. Path to VCF file of heterozygous SNPs.
#' @param p_onsv Character. Path to VCF file of heterozygous SNPs on
#'   SV-supporting reads.
#' @param p_sv Character. Path to VCF file of structural variants.
#' @param p_cnv Character. Path to CNV file (tab-delimited, FACETS format).
#' @param samp Character. Sample name written into the output table.
#' @param exper Character. Experiment identifier written into the output table.
#' @param chr_lst Character vector or NULL. Restrict analysis to these
#'   chromosomes (e.g. \code{c("chr1","chr2")}). Default \code{NULL} (all
#'   chromosomes).
#' @param flank_del Numeric. Maximum genomic distance (bp) for a deletion to be
#'   considered overlapping a translocation. Default \code{50}.
#' @param QUAL_tresh Numeric. Minimum SV quality score. Default \code{100}.
#' @param min_alt Numeric. Minimum SV-supporting read count. Default \code{2}.
#' @param tum_only Logical. Whether SVs were called from a tumor-only BAM.
#'   Default \code{FALSE}.
#' @param thresh Numeric (0, 1). Threshold for deciding SV/CNV event ordering
#'   in SVCF calculation. Default \code{0.1}.
#' @param flank_snp Numeric. Flanking window (bp) used when mapping SNPs to SVs.
#'   Default \code{500}.
#' @param flank_cnv Numeric. Flanking window (bp) used when assigning background
#'   CNV segments. Default \code{1000}.
#' @param run_clustering Logical. Run DP-GMM clustering on paired samples.
#'   Default \code{FALSE}.
#' @param pair_path Character. Path to a tab-delimited file with paired sample
#'   IDs (columns: \code{pre_BAT}, \code{on_BAT}). Required when
#'   \code{run_clustering = TRUE}.
#' @param pur_path Character. Path to a tab-delimited purity file (must contain
#'   columns \code{sample} and \code{purity}). Required when
#'   \code{run_clustering = TRUE}.
#' @param data_dir Character. Root directory containing SVCFit output BED files.
#'   Expected layout: \code{<data_dir>/COMBAT/SVCFit_output/<sample_ID>.bed}.
#'   Required when \code{run_clustering = TRUE}.
#' @param pair_num Integer. Which pair index to extract final clone CCFs for
#'   (used as input to tree building). Default \code{1L}.
#' @param exclude_pairs Integer vector. Pair indices to exclude from clustering.
#'   Default \code{integer(0)} (no pairs excluded).
#' @param Kmax Integer. Maximum number of DP-GMM components. Default \code{10L}.
#' @param n_steps Integer. Number of warm-start EM steps for convergence
#'   diagnostics. Default \code{100L}.
#' @param thr_min_w Numeric. Minimum component weight threshold below which a
#'   component is considered inactive. Default \code{0.01}.
#' @param random_state Integer. Python random seed for reproducibility.
#'   Default \code{0L}.
#' @param concentration Numeric. Dirichlet process concentration parameter.
#'   Default \code{1}.
#' @param min_n Integer. Minimum number of events per cluster; smaller clusters
#'   are merged into the nearest larger one. Default \code{5L}.
#' @param run_tree Logical. Build a phylogenetic tree from cluster CCFs.
#'   Requires \code{run_clustering = TRUE}. Default \code{FALSE}.
#' @param lineage_precedence_thresh Numeric. Maximum allowed CCF ratio for a
#'   child-to-parent edge to be kept in the tree graph. Default \code{0.2}.
#' @param sum_filter_thresh Numeric. Maximum allowed excess of children CCF
#'   sum over parent CCF. Default \code{0.2}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{svcf}}{data.frame. Full SVCF result table from
#'     \code{\link{calc_svcf}}.}
#'   \item{\code{cluster}}{(only when \code{run_clustering = TRUE}) List
#'     returned by \code{\link{cluster_data}}: \code{[[1]]} cluster assignments,
#'     \code{[[2]]} SV-to-cluster mapping, \code{[[3]]} clone CCF table.}
#'   \item{\code{tree}}{(only when \code{run_tree = TRUE}) List returned by
#'     \code{\link{build_tree}}: \code{[[1]]} best spanning tree edges,
#'     \code{[[2]]} CCF matrix, \code{[[3]]} igraph tree plot.}
#' }
#'
#' @examples
#' \dontrun{
#' # SVCF inference only
#' res <- run_svcfit(
#'   p_het  = "het_snp.vcf",
#'   p_onsv = "het_on_sv.vcf",
#'   p_sv   = "sv.vcf",
#'   p_cnv  = "cnv.txt",
#'   samp   = "SAMPLE_01",
#'   exper  = "EXP_01"
#' )
#' head(res$svcf)
#'
#' # SVCF inference + clustering + tree building
#' res <- run_svcfit(
#'   p_het          = "het_snp.vcf",
#'   p_onsv         = "het_on_sv.vcf",
#'   p_sv           = "sv.vcf",
#'   p_cnv          = "cnv.txt",
#'   samp           = "SAMPLE_01",
#'   exper          = "EXP_01",
#'   run_clustering = TRUE,
#'   pair_path      = "pairs.txt",
#'   pur_path       = "purity.txt",
#'   pair_num       = 1L,
#'   run_tree       = TRUE
#' )
#' res$tree[[3]]  # igraph plot of best tree
#' }
#'
#' @export
run_svcfit <- function(
    # --- SVCF inference (required) -------------------------------------------
    p_het,
    p_onsv,
    p_sv,
    p_cnv,
    samp,
    exper,

    # --- SVCF inference options -----------------------------------------------
    chr_lst    = NULL,
    flank_del  = 50,
    QUAL_tresh = 100,
    min_alt    = 2,
    tum_only   = FALSE,
    thresh     = 0.1,
    flank_snp  = 500,
    flank_cnv  = 1000,

    # --- Clustering (optional) ------------------------------------------------
    run_clustering = FALSE,
    pair_path      = NULL,
    pur_path       = NULL,
    data_dir       = NULL,
    pair_num       = 1L,
    exclude_pairs  = integer(0),
    Kmax           = 10L,
    n_steps        = 100L,
    thr_min_w      = 0.01,
    random_state   = 0L,
    concentration  = 1,
    min_n          = 5L,

    # --- Tree building (optional) ---------------------------------------------
    run_tree                  = FALSE,
    lineage_precedence_thresh = 0.2,
    sum_filter_thresh         = 0.2
) {

  # ---- Input validation ------------------------------------------------------
  if (run_tree && !run_clustering) {
    stop("run_tree = TRUE requires run_clustering = TRUE.")
  }
  if (run_clustering) {
    if (is.null(pair_path)) {
      stop("pair_path must be provided when run_clustering = TRUE.")
    }
    if (is.null(pur_path)) {
      stop("pur_path must be provided when run_clustering = TRUE.")
    }
    if (is.null(data_dir)) {
      stop("data_dir must be provided when run_clustering = TRUE.")
    }
  }

  # ---- Phase 1: Extract data -------------------------------------------------
  message("[SVCFit] Extracting data and parsing VCFs ...")
  extracted <- extract_info(
    p_het      = p_het,
    p_onsv     = p_onsv,
    p_sv       = p_sv,
    p_cnv      = p_cnv,
    chr_lst    = chr_lst,
    flank_del  = flank_del,
    QUAL_tresh = QUAL_tresh,
    min_alt    = min_alt,
    tum_only   = tum_only
  )

  data     <- extracted[[1]]   # list(het_snp, het_on_sv, sv, cnv)
  sv_info  <- extracted[[2]]   # data.frame from parse_sv_info
  sv_phase <- extracted[[4]]   # data.frame from parse_snp_on_sv

  # ---- Phase 2: Characterize SVs ---------------------------------------------
  message("[SVCFit] Characterizing SVs (CNV assignment, phasing) ...")
  anno_sv_cnv <- characterize_sv(
    sv_phase  = sv_phase,
    sv_info   = sv_info,
    cnv       = data$cnv,
    flank_snp = flank_snp,
    flank_cnv = flank_cnv
  )

  # ---- Phase 3: Calculate SVCF -----------------------------------------------
  message("[SVCFit] Calculating SVCF ...")
  svcf <- calc_svcf(
    anno_sv_cnv = anno_sv_cnv,
    sv_info     = sv_info,
    thresh      = thresh,
    samp        = samp,
    exper       = exper
  )

  result <- list(svcf = svcf)
  message("[SVCFit] SVCF inference complete. ", nrow(svcf), " SVs processed.")

  # ---- Phase 4 (optional): Clustering ----------------------------------------
  if (run_clustering) {
    message("[SVCFit] Running DP-GMM clustering on paired samples ...")
    cluster_out <- cluster_data(
      pair_path     = pair_path,
      pur_path      = pur_path,
      data_dir      = data_dir,
      Kmax          = Kmax,
      n_steps       = n_steps,
      thr_min_w     = thr_min_w,
      random_state  = random_state,
      concentration = concentration,
      min_n         = min_n,
      pair_num      = pair_num,
      exclude_pairs = exclude_pairs
    )
    result$cluster <- cluster_out
    message("[SVCFit] Clustering complete.")

    # ---- Phase 5 (optional): Tree building ------------------------------------
    if (run_tree) {
      message("[SVCFit] Building phylogenetic tree ...")
      clones <- cluster_out[[3]]  # clone CCF table for the selected pair
      if (nrow(clones) < 2) {
        warning("[SVCFit] Fewer than 2 clones found for pair_num = ", pair_num,
                "; skipping tree building.")
      } else {
        tree_out <- build_tree(
          clones                    = clones,
          lineage_precedence_thresh = lineage_precedence_thresh,
          sum_filter_thresh         = sum_filter_thresh
        )
        result$tree <- tree_out
        message("[SVCFit] Tree building complete.")
      }
    }
  }

  message("[SVCFit] Pipeline finished.")
  return(result)
}
