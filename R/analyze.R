#' main wrapper to run all functions
#'
#' @param sampan object of class 'character'. This object describes the sample name.
#' @param exper an object of class 'character'. This object describes the experiment number.
#' @param thresh an object of class 'numeric'. This object is strictly between 0 and 1, and describes the
#' threshold for SVCFit to decide whether an SV occurs before CNV or after.
#'
#' @returns sth
#' @export
#'
analyze <- function(samp, exper, thresh = 0.1, chr_lst=NULL) {
  print(c(samp,exper))
  print('now loading data')
  data      <- load_data(samp, exper, chr=chr_lst)
  print('now processing bnd')
  bnd_tmp   <- preproc_bnd(data$sv)
  del       <- process_del(data$sv, bnd_tmp, flank_del=50)
  bnd       <- annotate_bnd(bnd_tmp, del)
  print('now pharsing sv info')
  sv_info   <- parse_sv_info(data$sv, bnd, del, QUAL_tresh=100, min_alt=2)
  print('now pharsing het snp')
  snp_df    <- parse_het_snps(data$het_snp)
  print('now pharsing snp on sv')
  sv_phase  <- parse_snp_on_sv(data$het_on_sv, snp_df)
  print('now assign sv id to snp')
  assign_id <- assign_svids(sv_phase, sv_info, flank=500)
  print('now summarizing sv zygosity and phasing')
  sv_sum    <- sum_sv_info(sv_phase, assign_id, sv_info)
  print('now assinging cnv to sv')
  sv_cnv    <- assign_cnv(sv_sum, data$cnv)
  print('now creating cnv phasing')
  anno_sv_cnv     <- annotate_cnv(sv_cnv)
  print('now load the truth')
  truth <- load_truth(exper)
  print('now calculating svcf')
  final <- calc_svcf(anno_sv_cnv, sv_info, thresh=0.1, samp, exper)
  print('now attach truth to final')
  final_truth <- attach_truth(final, truth)
  print('now analyzing svclone')
  svc <- analyze_svclone(samp,exper, final_truth, sv_info)
  common_dat <- final_truth %>%
    filter(ID %in% svc$ID)
  return(list(final_truth, svc, common_dat))
}
