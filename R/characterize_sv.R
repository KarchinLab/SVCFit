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
  assign_id <- assign_svids(sv_phase, sv_info, flank=500)
  sv_sum    <- sum_sv_info(sv_phase, assign_id, sv_info)
  sv_cnv    <- assign_cnv(sv_sum, cnv)
  anno_sv_cnv     <- annotate_cnv(sv_cnv)
  sv_background <- assign_background_cnv(sv_sum, cnv, flank = 1000)
  sv_all <- left_join(anno_sv_cnv, sv_background)
  return(sv_all)
}
