#' Annotate the type and phasing of a CNV and calculate allele specific copy number for it.
#'
#' @param sv_cnv an object of class 'data frame'. This object stores the output from `assign_cnv`
#'
#' @returns sth
#' @export
#'
annotate_cnv <- function(sv_cnv) {
  anno_sv_cnv = sv_cnv %>%
    mutate(
      vaf = snp_alt/dep,
      tascn = round(vaf/(1-vaf),2),
      cn_type = case_when(
        cna > 2           ~ 'DUP',
        cna == 2 & minor==1 ~ 'norm',
        minor == 0        ~ 'DEL'
      ),
      cnv_phase = case_when(
        cn_type=='DUP' & tascn > 1  ~ 'pat',
        cn_type=='DUP' & tascn < 1  ~ 'mat',
        cn_type=='DEL' & tascn > 1   ~ 'mat',
        cn_type=='DEL' & tascn < 1   ~ 'pat',
        TRUE                         ~ 'pat'
      ),
      ASCN = case_when(
        cn_type=='DUP' & cnv_phase=='pat' ~ tascn,
        cn_type=='DUP' & cnv_phase=='mat' ~ 1/tascn,
        cn_type=='DEL' & cnv_phase=='pat' ~ tascn,
        cn_type=='DEL' & cnv_phase=='mat' ~ 1/tascn,
        TRUE                              ~ 1
      )
    )%>%
    arrange(POS)%>%
    select(CHROM, POS, ID, zygosity, sv_phase, cnv_phase, cncf, major, minor, cna, ASCN,cn_type, tascn, mate)
  return(anno_sv_cnv)
}
