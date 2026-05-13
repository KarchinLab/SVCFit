#' @import dplyr
#' @import tidyr
#' @importFrom ggplot2 ggplot aes geom_bar geom_line geom_point geom_hline
#'   coord_polar scale_fill_manual theme_void theme element_text facet_wrap
#'   labs theme_classic
#' @importFrom stringr str_detect
#' @importFrom GenomicRanges GRanges findOverlaps pintersect width seqnames mcols
#'   start end
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom grDevices colorRampPalette colors
#' @importFrom graphics par
#' @importFrom utils read.table read.delim
#' @importFrom stats dist setNames
#' @importFrom reticulate import py_config
"_PACKAGE"

# Silence R CMD check NOTEs about non-standard evaluation in dplyr/tidyr
# pipelines (bare column names referenced as if they were variables).
utils::globalVariables(c(
  ".", "AD", "ALT", "ASCN", "CHROM", "DEP", "END", "END_cluster", "END_raw",
  "FILTER", "FORMAT", "ID", "INFO", "POS", "POS_cluster", "POS_raw",
  "Proportion", "QUAL", "REF", "Subclone", "a_count", "allele",
  "best_sv_alt", "best_sv_ref", "bi_directional", "bkg_cnv", "bnd_group",
  "bulk", "c_birth", "c_die", "c_last_seen", "ccf", "ccube_ccf1",
  "ccube_ccf2", "child", "chr2", "chrom", "classification", "clone",
  "closest_keeper", "cluster", "cluster_merged", "cluster_num", "cn_type",
  "cna", "cncf", "cnv_phase", "connected", "delta", "dep", "donor",
  "doubleBreakPtsRes", "edge", "event_id", "f_day85", "f_pre", "final_svcf",
  "grp_left", "grp_right", "info", "iter", "key", "lcn.em", "left", "len",
  "len1", "len2", "len_dif", "lower_bound", "mPOS", "major", "match_chrom",
  "match_pos", "match_pos_id", "mate", "minor", "mutation_id", "nPOS",
  "n_in_clust", "ncluster", "new_on", "new_pre", "no_snp", "num", "on_BAT",
  "onsv_alt", "onsv_alt2", "onsv_ref", "overlap", "p_birth", "p_die",
  "p_last_seen", "pair", "parent", "pos2", "post_c", "post_ccf",
  "post_center", "pre_BAT", "pre_c", "pre_ccf", "pre_center", "prob",
  "purity", "r_2", "r_bar", "raw_svcf", "receiver", "reverse_edge",
  "reversed_connected", "right", "s1", "s2", "samp", "sample_ID",
  "se_day85", "se_pre", "second", "shared", "sid", "side", "snp_alt",
  "snp_ref", "ss1", "ss2", "stage", "sv_alt", "sv_len", "sv_phase",
  "sv_ref", "sv_zy", "tascn", "tmpid", "tumor", "type", "use_mate",
  "v_sorted", "vaf", "var_day85", "var_pre", "vio_birth", "w_day85",
  "w_pre", "z1", "z2", "zygosity", "zygosity_snp", "zygosity_tmp"
))
