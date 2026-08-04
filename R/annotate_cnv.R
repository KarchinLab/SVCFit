#' Annotate the type and phasing of a CNV and calculate the allele copy ratio.
#'
#' Ploidy-aware version. With \code{hemizygous_chr = NULL} the classification reduces exactly to the
#' original diploid tests. See chrX_rerun/MATH.md, change 3.
#'
#' On a hemizygous chromosome there are no heterozygous germline SNPs, so \code{tascn} and
#' \code{ASCN} (Eq. 6) are undefined. They are left as NA and the downstream correction is skipped
#' rather than being fed a fabricated value.
#'
#' @param sv_cnv data.frame. Output of `assign_cnv`.
#' @param hemizygous_chr character vector or NULL. Chromosomes that are single-copy in this
#'   subject's germline.
#'
#' @return data.frame with the original columns plus \code{pl}.
#' @export
#'
annotate_cnv <- function(sv_cnv, hemizygous_chr = NULL) {
  anno_sv_cnv <- sv_cnv %>%
    mutate(
      pl    = local_ploidy(CHROM, hemizygous_chr),
      vaf   = snp_alt / dep,
      tascn = round(vaf / (1 - vaf), 2),

      ## no heterozygous SNPs exist on a hemizygous chromosome, so the allele copy ratio is not
      ## estimable there; do not let a spurious value through
      tascn = ifelse(pl == 1L, NA_real_, tascn),

      ## ploidy-aware classification. For pl == 2 these are the original tests.
      cn_type = classify_cn(cna, minor, pl),

      cnv_phase = case_when(
        pl == 1L                     ~ NA_character_,
        cn_type == 'DUP' & tascn > 1 ~ 'pat',
        cn_type == 'DUP' & tascn < 1 ~ 'mat',
        cn_type == 'DEL' & tascn > 1 ~ 'mat',
        cn_type == 'DEL' & tascn < 1 ~ 'pat',
        TRUE                         ~ 'pat'
      ),

      ASCN = case_when(
        pl == 1L                                ~ NA_real_,
        cn_type == 'DUP' & cnv_phase == 'pat'   ~ tascn,
        cn_type == 'DUP' & cnv_phase == 'mat'   ~ 1 / tascn,
        cn_type == 'DEL' & cnv_phase == 'pat'   ~ tascn,
        cn_type == 'DEL' & cnv_phase == 'mat'   ~ 1 / tascn,
        TRUE                                    ~ 1
      )
    ) %>%
    arrange(POS) %>%
    select(CHROM, POS, ID, zygosity, sv_phase, cnv_phase, cncf, major, minor, cna,
           ASCN, cn_type, tascn, no_snp, mate, pl)

  return(anno_sv_cnv)
}
