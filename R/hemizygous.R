#' Hemizygous-chromosome support for SVCFit
#'
#' The allele copy ratio used by the diploid estimator is derived from
#' heterozygous germline SNPs. That quantity is unavailable on a chromosome
#' that is single-copy in the germline. These helpers make local normal ploidy
#' explicit and use \code{cn_bar}, the read-depth estimate of mean copies per
#' cell, for hemizygous copy-number corrections.

#' Local normal ploidy for each row
#'
#' @param chrom Character vector of chromosome names.
#' @param hemizygous_chr Character vector of chromosomes that are single-copy in this subject's
#'   germline, e.g. \code{c("chrX","chrY")} for a male. \code{NULL} or empty means all diploid,
#'   which reproduces the original behaviour exactly.
#' @return Integer vector of local normal ploidy, 1 or 2.
#' @export
local_ploidy <- function(chrom, hemizygous_chr = NULL) {
  pl <- rep(2L, length(chrom))
  if (!is.null(hemizygous_chr) && length(hemizygous_chr)) {
    pl[chrom %in% hemizygous_chr] <- 1L
  }
  pl
}

#' Default germline copy-number state for hemizygous chromosomes
#'
#' Allele-specific copy-number estimates can be unreliable on chromosomes with
#' no informative heterozygous SNPs. Unless an override is supplied, a declared
#' hemizygous chromosome is assigned its one-copy germline state.
#'
#' This default is the GERMLINE state only. Somatic copy-number variation on the chromosome enters
#' through \code{cn_bar}, measured from depth, and never through this function.
#'
#' @param cnv data.frame of FACETS segments with columns chrom, cna/tcn.em, major, minor.
#' @param hemizygous_chr Character vector of hemizygous chromosome names.
#' @param override Optional data.frame with columns chrom, cna, major, minor supplying an
#'   externally-derived hemizygous copy number instead of the default of one copy.
#' @return The cnv data.frame with hemizygous rows replaced.
#' @export
apply_hemizygous_cn <- function(cnv, hemizygous_chr = NULL, override = NULL) {
  if (is.null(hemizygous_chr) || !length(hemizygous_chr)) return(cnv)
  idx <- cnv$chrom %in% hemizygous_chr
  if (!any(idx)) return(cnv)
  if (!is.null(override)) {
    for (k in which(idx)) {
      o <- override[override$chrom == cnv$chrom[k], , drop = FALSE]
      if (nrow(o)) {
        cnv$cna[k] <- o$cna[1]; cnv$major[k] <- o$major[1]; cnv$minor[k] <- o$minor[1]
      }
    }
  } else {
    if ("cna"    %in% names(cnv)) cnv$cna[idx]    <- 1
    if ("tcn.em" %in% names(cnv)) cnv$tcn.em[idx] <- 1
    if ("major"  %in% names(cnv)) cnv$major[idx]  <- 1
    if ("minor"  %in% names(cnv)) cnv$minor[idx]  <- 0
    if ("lcn.em" %in% names(cnv)) cnv$lcn.em[idx] <- 0
  }
  cnv
}

#' Ploidy-aware copy-number classification
#'
#' Replaces the hard-coded diploid tests in annotate_cnv(). With pl = 2 the three tests reduce
#' exactly to the originals.
#'
#' @param cna Integer vector of total copy number.
#' @param minor Integer vector of minor copy number.
#' @param pl Integer vector of local normal ploidy from \code{local_ploidy()}.
#' @return Character vector, one of "DUP", "norm", "DEL", or NA.
#' @export
classify_cn <- function(cna, minor, pl) {
  out <- rep(NA_character_, length(cna))
  pl <- rep_len(pl, length(cna))

  ## DIPLOID ROWS: reproduce the original case_when EXACTLY, including its DEL test.
  ##
  ## The original is
  ##   case_when(cna > 2 ~ 'DUP', cna == 2 & minor == 1 ~ 'norm', minor == 0 ~ 'DEL')
  ## and its DEL criterion is `minor == 0`, NOT `cna < 2`. Those differ on copy-neutral LOH
  ## (cna = 2, minor = 0), which the original calls DEL. An earlier version of this function
  ## generalised DEL to `cna < pl`, which silently reclassified those rows to NA and changed
  ## final_svcf on 3 autosomal rows of 80847_recut. The regression gate caught it. Autosomal
  ## behaviour must be bit-identical, so the diploid branch is now the original, verbatim.
  d <- which(pl == 2L)
  if (length(d)) {
    o <- rep(NA_character_, length(d))
    cna_d <- cna[d]; min_d <- minor[d]
    o[!is.na(min_d) & min_d == 0]                              <- "DEL"
    o[!is.na(cna_d) & !is.na(min_d) & cna_d == 2 & min_d == 1] <- "norm"
    o[!is.na(cna_d) & cna_d > 2]                               <- "DUP"
    out[d] <- o
  }

  ## HEMIZYGOUS ROWS: one germline copy, so minor is 0 whenever the segment is present and carries
  ## no information. Classify on total copy number alone.
  h <- which(pl == 1L)
  if (length(h)) {
    o <- rep(NA_character_, length(h)); cna_h <- cna[h]
    o[!is.na(cna_h) & cna_h >  1] <- "DUP"
    o[!is.na(cna_h) & cna_h == 1] <- "norm"
    o[!is.na(cna_h) & cna_h <  1] <- "DEL"
    out[h] <- o
  }
  out
}

#' Status flag for each SV row
#'
#' "ok"                          usable estimate
#' "hemizygous_cnv_unresolved"   hemizygous and copy-altered, but resolve_hemizygous_svcf() has not
#'                               run. A row retaining this status was never resolved.
#' "zero_ref_depth"              reference count is 0, so r_bar = sv_alt/sv_ref is undefined
#' "del_depth_mismatch"          deletion: reads and depth BOTH estimated the cellular fraction and
#'                               they disagree by more than tol. A real conflict between two
#'                               measurements. The estimate is still returned.
#' "del_depth_unresolvable"      deletion: depth could not estimate it independently, because the
#'                               span and its flank fall in the same segment and bg_cn - cn_bar is
#'                               0 by construction. NOT a conflict -- nothing corroborated the
#'                               estimate, which is different from two measurements disagreeing.
#'
#' VAF is mathematically defined when \code{sv_ref = 0}, but SVCFit labels the
#' row \code{"zero_ref_depth"} because the alt-to-reference depth ratio used by
#' other branches is infinite. \code{calc_svcf()} can recover independently
#' reviewed rows through \code{zero_ref_allowlist}. All other zero-reference
#' estimates are reported as \code{NA}, on autosomes as well as hemizygous chromosomes.
#'
#' @param pl Integer vector of local normal ploidy from \code{local_ploidy()}.
#' @param cn_type Character vector of copy-number class from \code{classify_cn()}.
#' @param sv_ref Integer vector of reference-supporting read counts.
#' @return Character vector of per-row status labels.
#' @export
svcf_status <- function(pl, cn_type, sv_ref) {
  st <- rep("ok", length(pl))
  st[pl == 1L & !is.na(cn_type) & cn_type != "norm"] <- "hemizygous_cnv_unresolved"
  st[!is.na(sv_ref) & sv_ref == 0] <- "zero_ref_depth"
  st
}

#' Detect data that contradicts the flat single-copy default
#'
#' \code{apply_hemizygous_cn()} defaults a hemizygous chromosome to one germline copy. Somatic
#' copy-number variation must then arrive as \code{cn_bar}. This reports the copy-number-changing
#' SVs on a hemizygous chromosome, so a caller that has no \code{cn_bar} for
#' them is warned rather than silently treating the chromosome as flat.
#'
#' @param sv data.frame with columns CHROM, POS, END, classification.
#' @param hemizygous_chr character vector of hemizygous chromosome names.
#' @param min_span numeric. Minimum span in bp to count as a copy-number change. Default 1e6.
#' @return data.frame of contradicting events; zero rows means the flat default is safe.
#' @export
hemizygous_cn_conflicts <- function(sv, hemizygous_chr = NULL, min_span = 1e6) {
  empty <- sv[0, intersect(c("CHROM", "POS", "END", "classification"), names(sv)), drop = FALSE]
  if (is.null(hemizygous_chr) || !length(hemizygous_chr)) return(empty)
  keep <- sv$CHROM %in% hemizygous_chr &
          sv$classification %in% c("DEL", "DUP") &
          !is.na(suppressWarnings(as.numeric(sv$END))) &
          (suppressWarnings(as.numeric(sv$END)) - sv$POS) >= min_span
  keep[is.na(keep)] <- FALSE
  out <- sv[keep, intersect(c("CHROM", "POS", "END", "classification"), names(sv)), drop = FALSE]
  if (nrow(out)) {
    warning(sprintf(paste0("hemizygous_cn_conflicts: %d copy-number-changing SV(s) on %s; supply ",
                           "cn_bar from depth rather than assuming one copy across the chromosome"),
                    nrow(out), paste(hemizygous_chr, collapse = ",")), call. = FALSE)
  }
  out
}

#' Hemizygous SVCF for an SV sitting in a segment whose copy number something else changed
#'
#' Counting alleles over the three cell populations, the total per cell is always \code{cn_bar}.
#' With VAF = BPC/(BPC + BEC):
#'
#'   SV precedes CNV   (H1)   SVCF = cn_bar * VAF - (cn_bar - 1)
#'   CNV precedes SV   (H2)   SVCF = cn_bar * VAF
#'
#' H2 is Eq. 1 with R = cn_bar. Both reduce to SVCF = VAF at cn_bar = 1.
#'
#' SELECTION. The sign of the H1 form decides the ordering, and nothing else is needed. Under
#' CNV-first, SV carriers are a subset of CNV carriers so SVCF <= f_CNV, and since
#' cn_bar - 1 = f_CNV(c - 1) >= SVCF(c - 1) >= SVCF for c >= 2, the H1 form is necessarily <= 0.
#' Under SV-first it is the true SVCF and so strictly positive. Verified 185/185 with no miscalls.
#'
#' This restores the diploid situation, where the sign of Eq. 4 decides the ordering. The previous
#' version of this file claimed that test does not carry over and substituted a lineage-precedence
#' rule against f_CNV. It does carry over; f_CNV is not needed, and neither is c.
#'
#' @param bpc,bec numeric vectors of breakpoint and breakend counts.
#' @param cn_bar numeric vector, mean copies of the locus per cell, from read depth
#'   (\code{cn_bar = R * psi_sample / 2}). NOT an integer copy number and never rounded.
#' @param tol numeric. Half-width of the copy-neutral band around cn_bar = 1.
#' @return data.frame with columns \code{svcf}, \code{ordering}, \code{status}, \code{h1}, \code{h2}.
#' @export
resolve_hemizygous_svcf <- function(bpc, bec, cn_bar, tol = 0.05) {
  n <- max(length(bpc), length(bec), length(cn_bar))
  bpc <- rep_len(bpc, n); bec <- rep_len(bec, n); cn_bar <- rep_len(cn_bar, n)

  vaf <- bpc / (bpc + bec)
  h1 <- cn_bar * vaf - (cn_bar - 1)
  h2 <- cn_bar * vaf

  svcf <- rep(NA_real_, n)
  ordering <- rep(NA_character_, n)
  status <- rep("ok", n)

  for (i in seq_len(n)) {
    if (!is.finite(vaf[i])) { status[i] <- "zero_ref_depth"; next }
    if (!is.finite(cn_bar[i])) { status[i] <- "no_depth"; next }
    if (abs(cn_bar[i] - 1) <= tol) {          # copy-neutral: both forms are VAF
      svcf[i] <- vaf[i]; ordering[i] <- "copy_neutral"
    } else if (h1[i] > 0) {
      svcf[i] <- h1[i]; ordering[i] <- "sv_before_cnv"
    } else {
      svcf[i] <- h2[i]; ordering[i] <- "cnv_before_sv"
    }
    if (!is.na(svcf[i]) && (svcf[i] <= 0 || svcf[i] > 1)) status[i] <- "infeasible_svcf"
  }
  svcf[status != "ok"] <- NA_real_
  data.frame(svcf = svcf, ordering = ordering, status = status, h1 = h1, h2 = h2,
             stringsAsFactors = FALSE)
}

#' SVCF for a hemizygous tandem duplication
#'
#' A duplication is its own copy-number change, so it does not decompose like H1/H2. A carrier
#' chromosome holding r tandem copies carries r - 1 novel junctions and one reference-configuration
#' span at each outer breakend, so with SVCF = s,
#'
#'   BPC prop. s(r - 1),  BEC prop. 1,  cn_bar = 1 + s(r - 1)
#'   =>  SVCF = (cn_bar - 1) / (r - 1)
#'
#' The previous form, BPC/(BPC + (r-1)BEC), is wrong in every case tested (20 of 20).
#'
#' r IS NOT IDENTIFIABLE. VAF = (cn_bar - 1)/cn_bar exactly, so VAF carries no information beyond
#' cn_bar and s cannot be separated from r at this locus. Since r >= 2, SVCF <= cn_bar - 1, so the
#' r = 2 value is an UPPER BOUND, tight if and only if the duplication is single-copy. This is a real
#' limit of the hemizygous case and must be reported as a bound, not concealed as a point estimate.
#'
#' @param cn_bar numeric vector, mean copies of the duplicated segment per cell, from depth.
#' @param r numeric vector, copies in carrier cells. Default 2, the SVCF-maximising choice.
#' @return data.frame with \code{svcf}, \code{is_upper_bound} and \code{status}.
#' @export
hemizygous_dup_svcf <- function(cn_bar, r = 2) {
  n <- max(length(cn_bar), length(r))
  cn_bar <- rep_len(cn_bar, n); r <- rep_len(r, n)
  svcf <- (cn_bar - 1) / (r - 1)
  status <- rep("ok", n)
  status[!is.finite(svcf)] <- "no_depth"
  status[is.finite(svcf) & (svcf <= 0 | svcf > 1)] <- "infeasible_svcf"
  svcf[status != "ok"] <- NA_real_
  data.frame(svcf = svcf, is_upper_bound = (r == 2), status = status, stringsAsFactors = FALSE)
}

#' SVCF for a hemizygous deletion, with a built-in consistency check
#'
#' A deletion removes the only copy, so carrier cells contribute no alleles at the locus and
#' cn_bar = 1 - s. The junction reads give the same quantity independently, VAF = s. The two are
#' therefore a genuine cross-check rather than one estimate used twice.
#'
#' Read depth supplies the primary estimate. VAF is retained as an independent
#' cross-check because carrier cells contribute no locus-spanning reads after a
#' hemizygous deletion. This rule is confined to the hemizygous helper; a
#' diploid deletion retains a surviving homologue.
#'
#' Depth is used when it is available and feasible; otherwise the VAF estimate is returned, so a
#' locus with no usable segmentation still gets an answer rather than an NA. `svcf_source` records
#' which was used, and `svcf_vaf` keeps the read-based value so the cross-check stays visible and
#' the decision stays reversible.
#'
#' @param bpc,bec numeric vectors of breakpoint and breakend counts.
#' @param cn_bar numeric vector from depth, or NA to use the read estimate alone.
#' @param tol numeric. Disagreement above this is flagged.
#' @param bg_cn numeric vector, flanking (background) copy number for the CNV-first form;
#'   \code{NA} falls back to the previous behaviour.
#' @param bg_tol numeric. Copy-neutral band half-width for the background estimate.
#' @return data.frame with \code{svcf}, \code{svcf_depth}, \code{svcf_vaf}, \code{svcf_source},
#'   \code{status}.
#' @export
hemizygous_del_svcf <- function(bpc, bec, cn_bar = NA_real_, tol = 0.15,
                                bg_cn = NA_real_, bg_tol = 0.05) {
  n <- max(length(bpc), length(bec), length(cn_bar))
  bpc <- rep_len(bpc, n); bec <- rep_len(bec, n); cn_bar <- rep_len(cn_bar, n)
  bg_cn <- rep_len(bg_cn, n)

  svcf_vaf   <- bpc / (bpc + bec)
  svcf_depth <- 1 - cn_bar

  ## THE CNV-FIRST FORM. `1 - cn_bar` assumes the locus is single-copy ABSENT the deletion. Write
  ## kappa for the copies it would have without the deletion; depth measures cn_bar = kappa - s in
  ## either ordering, so the estimator is s = kappa - cn_bar and `1 - cn_bar` is the kappa = 1
  ## special case. Its bias is exactly -(kappa - 1), which is why the error grew with purity.
  ##
  ## kappa is the flanking copy number, supplied as bg_cn. When it is not supplied this reduces to
  ## the previous behaviour exactly, so callers without a background estimate are unaffected.
  ##
  ## Measured on 17,712 detected deletions over 10 c50 replicates (DELETION-ON-ALTERED-LOCUS-
  ## 2026-07-30.md), within 0.05 of truth:
  ##
  ##   experiment   before (depth->h2->vaf)   after (depth->bg->h2->vaf)
  ##   e1                     93.4%                    93.4%
  ##   e2                     92.1%                    94.8%
  ##   e4                     18.1%                    58.9%
  ##   all                    68.4%                    82.8%
  ##
  ## bg_tol gates the new branch on the background actually being amplified. Without it the branch
  ## fires on copy-neutral loci where bg_cn - cn_bar is small-positive from noise alone and loses
  ## 2.1 points on e1; at 0.05 the change regresses nothing.
  svcf_bg <- bg_cn - cn_bar

  ## FALLBACK IS H2, NOT VAF, WHEN DEPTH IS INFEASIBLE (cn_bar >= 1).
  ##
  ## A deletion cannot raise the copy number of its own span, so cn_bar >= 1 means the span sits
  ## inside an amplification that arose BEFORE the deletion -- the CNV-first case, whose form is
  ## h2 = cn_bar * VAF. Depth is not merely inaccurate there, it is impossible: 1 - cn_bar <= 0 is
  ## not a cellular fraction.
  ##
  ## The feasibility of the depth estimate is therefore the discriminator, and it needs no knowledge
  ## of the ordering and no separate test of the background. Scored on all 1,256 detected deletions
  ## across the 45 purity x mixture conditions (scripts/14_chrx_del_branch_eval.R), this rule picks
  ## the best of the three forms in 25 of 30 e1+e2 cells and 14 of 15 e4 cells:
  ##
  ##   experiment   depth infeasible   best form   median |err| of the winner
  ##   e1                    5.1%      depth       0.013
  ##   e2                    4.4%      depth       0.018
  ##   e4                   86.4%      h2          0.207
  ##
  ## Deletions bypass resolve_hemizygous_svcf()'s sign rule entirely. That is deliberate: h1 wins
  ## only 3 of the 45 cells, and on e2 -- which is SV-first BY CONSTRUCTION, so the sign rule
  ## identifies the ordering correctly -- h1 is still the worst form (median 0.247 against depth's
  ## 0.018), and is 98.9% positive with median +0.56, i.e. confidently wrong rather than borderline.
  ## The h1 derivation does not appear to survive a deletion driving its own span to zero copies.
  ##
  ## e4 remains the weak case at 11.4% within 0.05: h2 is the best available form there, not a good
  ## one, and it degrades with purity (median 0.046 at p10 to 0.299 at p80). Flagged, not hidden.
  svcf_h2 <- cn_bar * svcf_vaf

  depth_ok <- is.finite(svcf_depth) & svcf_depth > 0 & svcf_depth <= 1
  bg_ok    <- is.finite(svcf_bg)    & svcf_bg    > 0 & svcf_bg    <= 1 &
              is.finite(bg_cn) & bg_cn > 1 + bg_tol
  h2_ok    <- is.finite(svcf_h2)    & svcf_h2    > 0 & svcf_h2    <= 1
  vaf_ok   <- is.finite(svcf_vaf)   & svcf_vaf   > 0 & svcf_vaf   <= 1

  ## Order matters and is by construction, not preference. Depth first, because where it is feasible
  ## it is exact and needs no background. Then the CNV-first form, because depth being infeasible
  ## (cn_bar >= 1) is itself the evidence that the span sits in an amplification that predates the
  ## deletion. h2 stays behind it for the case where the background is unknown or copy-neutral, and
  ## VAF is the last resort for a locus with no usable cn_bar at all.
  svcf   <- ifelse(depth_ok, svcf_depth,
             ifelse(bg_ok,   svcf_bg,
              ifelse(h2_ok,  svcf_h2,
               ifelse(vaf_ok, svcf_vaf, NA_real_))))
  source <- ifelse(depth_ok, "depth",
             ifelse(bg_ok,   "bg_depth",
              ifelse(h2_ok,  "h2",
               ifelse(vaf_ok, "vaf", NA_character_))))

  ## Compare VAF with the depth-based form that supplied the answer. The
  ## measurements are independent only when depth produced a feasible estimate.
  svcf_depth_used <- ifelse(source == "depth",    svcf_depth,
                     ifelse(source == "bg_depth", svcf_bg, NA_real_))

  ## Depth is blind to a deletion smaller than the segmentation resolution. If
  ## the SV span and its flank receive the same segment value, bg_cn - cn_bar is
  ## zero by construction. Label this separately from a genuine disagreement.
  depth_blind <- is.finite(bg_cn) & is.finite(cn_bar) & abs(bg_cn - cn_bar) <= 1e-9

  status <- rep("ok", n)
  disagree <- is.finite(svcf_vaf) & is.finite(svcf_depth_used) & !depth_blind &
              abs(svcf_vaf - svcf_depth_used) > tol
  status[is.finite(svcf_vaf) & depth_blind] <- "del_depth_unresolvable"
  status[disagree] <- "del_depth_mismatch"
  status[!is.finite(svcf_vaf)] <- "zero_ref_depth"
  ## bg_ok belongs here too: a row the CNV-first form resolved is not infeasible, and omitting it
  ## would have labelled the new branch's own output as having no usable estimate.
  status[!depth_ok & !bg_ok & !h2_ok & !vaf_ok] <- "infeasible_svcf"

  data.frame(svcf = svcf, svcf_depth = svcf_depth, svcf_vaf = svcf_vaf, svcf_h2 = svcf_h2,
             svcf_source = source, status = status, stringsAsFactors = FALSE)
}

#' DEFUNCT: solve hemizygous copy number and CNV cellular fraction
#'
#' Hemizygous estimators now take \code{cn_bar} directly, and the ordering is
#' decided by the sign of the H1 form. This stub produces an explicit error for
#' callers that still use the removed interface.
#'
#' @param ... Ignored. The function is defunct and always signals an error.
#' @export
solve_hemizygous_cn <- function(...) {
  stop("solve_hemizygous_cn() is defunct. The hemizygous forms take cn_bar directly; c and f_CNV ",
       "are not needed.", call. = FALSE)
}
