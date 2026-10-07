#' Annotate the type and phasing of a CNV and calculate the allele copy ratio.
#'
#' Ploidy-aware version. With \code{hemizygous_chr = NULL} the classification
#' reduces to the diploid tests.
#'
#' On a hemizygous chromosome there are no heterozygous germline SNPs, so \code{acr_raw} and
#' \code{acr} (Eq. 6) are undefined. They are left as NA and the downstream correction is skipped
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
      acr_raw = round(vaf / (1 - vaf), 2),

      ## no heterozygous SNPs exist on a hemizygous chromosome, so the allele copy ratio is not
      ## estimable there; do not let a spurious value through
      acr_raw = ifelse(pl == 1L, NA_real_, acr_raw),

      ## ploidy-aware classification. For pl == 2 these are the original tests.
      cn_type = classify_cn(cna, minor, pl),

      cnv_phase = case_when(
        pl == 1L                     ~ NA_character_,
        cn_type == 'DUP' & acr_raw > 1 ~ 'pat',
        cn_type == 'DUP' & acr_raw < 1 ~ 'mat',
        cn_type == 'DEL' & acr_raw > 1 ~ 'mat',
        cn_type == 'DEL' & acr_raw < 1 ~ 'pat',
        TRUE                         ~ 'pat'
      ),

      acr = case_when(
        pl == 1L                                ~ NA_real_,
        cn_type == 'DUP' & cnv_phase == 'pat'   ~ acr_raw,
        cn_type == 'DUP' & cnv_phase == 'mat'   ~ 1 / acr_raw,
        cn_type == 'DEL' & cnv_phase == 'pat'   ~ acr_raw,
        cn_type == 'DEL' & cnv_phase == 'mat'   ~ 1 / acr_raw,
        TRUE                                    ~ 1
      )
    ) %>%
    arrange(POS) %>%
    select(CHROM, POS, ID, zygosity, sv_phase, cnv_phase, cncf, major, minor, cna,
           acr, cn_type, acr_raw, no_snp, mate, pl)

  return(anno_sv_cnv)
}
