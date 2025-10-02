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
analyze <- function(sample, experiment, threshold = 0.1, chr_list=NULL) {
  print(c(sample,experiment))
  print('now loading data')
  data      <- load_data(sample, experiment)
  print('now processing bnd')
  bnd_tmp   <- preproc_bnd(data$sv)
  del       <- process_del(data$sv, bnd_tmp)
  bnd       <- annotate_bnd(bnd_tmp, del)
  print('now pharsing sv info')
  sv_info   <- parse_sv_info(data$sv, bnd, del)
  print('now pharsing het snp')
  snp_df    <- parse_het_snps(data$het_snp)
  print('now pharsing snp on sv')
  sv_phase  <- parse_snp_on_sv(data$het_on_sv, snp_df)
  print('now assign sv id to snp')
  assign_id <- assign_svids(sv_phase, sv_info)
  print('now summarizing sv zygosity and phasing')
  sv_sum    <- sum_sv_info(sv_phase, assign_id, sv_info)
  print('now assinging cnv to sv')
  sv_cnv    <- assign_cnv(sv_sum, data$cnv)
  print('now creating cnv phasing')
  anno_sv_cnv     <- annotate_cnv(sv_cnv)
  print('now load the truth')
  truth <- load_truth(experiment)
  print('now calculating svcf')
  final <- calc_svcf(anno_sv_cnv, sv_info, thresh=0.1, sample, experiment)
  print('now attach truth to final')
  final_truth <- attach_truth(final, truth)
  print('now analyzing svclone')
  svc <- analyze_svclone(sample,experiment, final_truth, sv_info)
  common_dat <- final_truth %>%
    filter(ID %in% svc$ID)
  print('now svcfit w svclone read')
  w_svc_read=svcfit_w_svc(svc, final_truth, thresh=0.1, sample, experiment)
  return(list(final_truth, svc, common_dat, w_svc_read))
}
