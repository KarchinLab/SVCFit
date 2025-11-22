#' Characterize SV
#'
#' @param sv_phase object of class 'dataframe'. This object stores the SV phasing
#'  and zygosity information obtained from *parse_snp_on_sv()*.
#' @param sv_info object of class 'dataframe'. This object stores SV genomic locations
#' and types obtained from *parse_sv_info()*. 
#' @param cnv object of class 'dataframe'. This object stores CNV information
#'  obtained from *load_data()*. 
#'
#' @returns sth
#' @export
#'

characterize_sv <- function(sv_phase, sv_info, cnv){
  print('now assign sv id to snp')
  assign_id <- assign_svids(sv_phase, sv_info, flank=500)
  print('now summarizing sv zygosity and phasing')
  sv_sum    <- sum_sv_info(sv_phase, assign_id, sv_info)
  print('now assinging cnv to sv')
  sv_cnv    <- assign_cnv(sv_sum, cnv)
  print('now creating cnv phasing')
  anno_sv_cnv     <- annotate_cnv(sv_cnv)
  return(anno_sv_cnv)
}