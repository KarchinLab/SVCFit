#' Annotate the type and phasing of a CNV and calculate allele specific copy number for it.
#'
#' @param sv_cnv an object of class 'data frame'. This object stores the output from `assign_cnv`
#'
#' @return A data.frame with columns \code{CHROM}, \code{POS}, \code{ID},
#'   \code{zygosity}, \code{sv_phase}, \code{cnv_phase}, \code{cncf},
#'   \code{major}, \code{minor}, \code{cna}, \code{ASCN} (allele-specific copy
#'   number used in SVCF calculation), \code{cn_type} (\code{"DUP"},
#'   \code{"norm"}, or \code{"DEL"}), \code{tascn}, \code{no_snp}, and
#'   \code{mate}.
#' @export
#'
annotate_cnv <- function(sv_cnv) {
  anno_sv_cnv = sv_cnv %>%
    # label CNA
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
    select(CHROM, POS, ID, zygosity, sv_phase, cnv_phase, cncf, major, minor, cna, ASCN,cn_type, tascn, no_snp, mate)
  return(anno_sv_cnv)
}
