#' Calculate SVCF
#'
#' Ploidy-aware version. With \code{hemizygous_chr = NULL} every line reduces to the original
#' diploid code and results are bit-identical; this is asserted by
#' chrX_rerun/scripts/04_regression_autosomes.R. See chrX_rerun/MATH.md, as corrected by
#' CORRECTION-c-vs-cnbar-2026-07-29.md.
#'
#' UPDATED 2026-07-29. The hemizygous branch previously took the CNV cellular fraction
#' (\code{hemi_cncf}) and the integer copy number, and selected the SV-CNV ordering by a
#' lineage-precedence rule. Both quantities are gone. The corrected forms take the measured mean
#' copy number \code{cn_bar} directly, and the ordering is decided by the sign of the SV-first form.
#'
#' @param anno_sv_cnv data.frame. Output of `annotate_cnv`.
#' @param sv_info data.frame. Output of `parse_sv_info`.
#' @param thresh numeric in (0, 1). Threshold for deciding whether the SV occurs before or after the
#'   CNV.
#' @param samp character. Sample name.
#' @param exper character. Experiment identifier.
#' @param hemizygous_chr character vector or NULL. Chromosomes that are single-copy in this
#'   subject's germline, e.g. c("chrX", "chrY") for a male subject. NULL reproduces the original
#'   diploid-only behaviour.
#' @param hemi_cn_bar numeric or NULL. Mean copies of the locus per cell for each row on a
#'   hemizygous chromosome, from read depth (\code{cn_bar = R * psi_sample / 2}). This is a
#'   continuous measurement and is never rounded to an integer copy number. NULL leaves hemizygous
#'   copy-altered rows unresolved rather than guessing, which is countable downstream via
#'   \code{svcf_status}.
#' @param zero_ref_allowlist data.frame or NULL. Rows on a hemizygous chromosome with sv_ref = 0
#'   that BAM evidence shows are genuine clonal losses rather than reference dropout, so that
#'   SVCF = VAF = 1 is the correct answer. Columns \code{sample}, \code{chrom}, \code{pos}, and
#'   \code{verdict}; only rows with \code{verdict == "recover"} are used.
#'
#'   GENERATED, NEVER HAND-WRITTEN. Produce it with scripts/13_zero_ref_allowlist.sh, which derives
#'   the verdict from matched-normal soft-clipping calibrated against control loci. The distinction
#'   cannot be made from the bed alone: both suppression guards fire on sv_ref == 0 and the matched
#'   normal is what separates a somatic clonal loss from a germline or mapping artifact.
#'
#' @param hemi_dup_r numeric. Copies in carrier cells for a hemizygous tandem duplication. Default 2.
#'   r is not identifiable from a single locus, and r = 2 maximises SVCF, so the duplication values
#'   are upper bounds; \code{svcf_is_bound} marks them.
#'
#' @return data.frame as before, plus \code{pl} (local normal ploidy), \code{svcf_status},
#'   \code{sv_cnv_order} and \code{svcf_is_bound}. Rows whose status is not "ok" are returned with
#'   \code{final_svcf = NA} rather than dropped, so exclusions are countable instead of silent.
#' @export
#'
calc_svcf <- function(anno_sv_cnv, sv_info, thresh = 0.1, samp, exper,
                      hemizygous_chr = NULL, hemi_cn_bar = NULL, hemi_dup_r = 2,
                      zero_ref_allowlist = NULL, hemi_bg_cn = NULL) {
  final <- inner_join(anno_sv_cnv, sv_info) %>%
    mutate(
      ## local normal ploidy: 2 on autosomes, 1 on a declared hemizygous chromosome
      pl = local_ploidy(CHROM, hemizygous_chr),

      ## Eq. 1 conversion applied to the read-depth ratio. Was 2*(alt+ref)/ref.
      r_bar = pl * (sv_alt + sv_ref) / sv_ref,

      ## extra copies above the local normal ploidy. Was (major+minor-2).
      r_2 = ifelse((major + minor - pl) < 1,
                   r_bar * sv_alt / (sv_alt + sv_ref),
                   (major + minor - pl)),

      ## Eq. 1. Was 2*alt/(alt+ref); on a hemizygous chromosome SVCF = VAF.
      raw_svcf = ifelse(classification == "DUP",
                        (r_bar * sv_alt / (sv_alt + sv_ref)) / r_2,
                        pl * sv_alt / (sv_alt + sv_ref)),

      ## the homozygous halving is a two-copy concept: an SV on both copies of a diploid locus.
      ## With one copy there is nothing to halve.
      raw_svcf = ifelse(zygosity == 'hom' & pl == 2L, raw_svcf / 2, raw_svcf),

      ## Eqs. 4 and 5. Both are built on ACR and are defined for pl == 2 only.
      s1  = (2 * sv_alt - sv_ref * (ASCN - 1)) / (sv_alt + sv_ref),
      s2  = ((sv_alt + sv_alt * ASCN) / (sv_alt + sv_ref)) %% 2,
      ss1 = ifelse(zygosity == "hom", 0.5 * s1, s1),
      ss2 = ifelse(zygosity == "hom", 0.5 * s2, s2),
      s1  = ifelse(pl == 1L, NA_real_, s1),
      s2  = ifelse(pl == 1L, NA_real_, s2),
      ss1 = ifelse(pl == 1L, NA_real_, ss1),
      ss2 = ifelse(pl == 1L, NA_real_, ss2),

      ## read-count based rule selection between Eqs. 4 and 5, unchanged for pl == 2
      final_svcf = ifelse(cn_type == "DEL" | ss1 <= thresh, ss2, ss1),
      final_svcf = ifelse(classification == "COPBND" & bkg_cnv == 'norm',
                          round(raw_svcf, 2), round(final_svcf, 2)),
      final_svcf = ifelse(bkg_cnv == 'norm' & classification == 'DEL',
                          round(raw_svcf, 2), round(final_svcf, 2)),
      final_svcf = ifelse(bkg_cnv == 'norm' & classification == 'DUP',
                          round(raw_svcf, 2), round(final_svcf, 2)),

      ## a copy-neutral hemizygous locus needs no CNV correction: SVCF = VAF (MATH.md change 1).
      ## Duplications are excluded here and handled below: see the next comment.
      final_svcf = ifelse(pl == 1L & !is.na(cn_type) & cn_type == 'norm' & classification != 'DUP',
                          round(raw_svcf, 2), final_svcf),

      ## A HEMIZYGOUS TANDEM DUPLICATION IS ITS OWN COPY-NUMBER CHANGE, so it is never genuinely
      ## copy-neutral at its own locus and the branch above must not claim it.
      ##
      ## The `raw_svcf` DUP formula also collapses here. With apply_hemizygous_cn() asserting
      ## major = 1, minor = 0 and pl = 1, the "extra copies above local ploidy" term
      ## (major + minor - pl) is 0, so r_2 falls back to r_bar*VAF and
      ##     raw_svcf = (r_bar*VAF) / (r_bar*VAF) = 1
      ## identically, for ANY input. Measured on the 2026-07-29 cohort run: all 13 chrX DUP rows
      ## returned exactly 1.00 while their VAFs ranged 0.224 to 0.506, which would have put 13
      ## spurious clonal duplications onto the trees. On autosomes the fallback is unreachable for
      ## a real duplication, because major + minor - 2 >= 1 there.
      ##
      ## The correct form needs cn_bar (H3, SVCF = (cn_bar - 1)/(r - 1)) and is applied in the
      ## block below when it is supplied. Without cn_bar the cellular fraction is NOT estimable,
      ## so mark it and leave it NA rather than emit a number.
      final_svcf  = ifelse(pl == 1L & classification == 'DUP', NA_real_, final_svcf),

      final_svcf = ifelse(is.na(final_svcf), round(raw_svcf, 2), final_svcf),

      ## Make every exclusion explicit and countable rather than a silent filter.
      ##
      ## SCOPE CONTAINMENT (RK, 2026-07-29). The status column is computed for every row, but the
      ## NA that follows from it is applied ONLY on a hemizygous chromosome. Autosomal output
      ## therefore stays bit-identical to the submitted run.
      ##
      ## Why this is not merely cosmetic: 203 autosomal rows across the cohort (3.9%) have
      ## sv_ref == 0, and the submitted pipeline reports final_svcf of 1 or 2 for them. Two is not
      ## a cellular fraction, so those values are wrong, and this patch would have corrected them
      ## to NA. That correction is real but out of scope for a revision about chrX: it widens the
      ## diff, forces the autosomal figures to be regenerated, and invites a reviewer question
      ## unrelated to the change under review. It is recorded in GATE-RESULT-2026-07-29.md and
      ## belongs in its own piece of work. `svcf_status` still labels these rows, so they remain
      ## countable without being altered.
      svcf_status = svcf_status(pl, cn_type, sv_ref),
      svcf_status = ifelse(pl == 1L & classification == 'DUP' & svcf_status == "ok",
                           "hemizygous_dup_needs_cn_bar", svcf_status),
      final_svcf  = ifelse(pl == 1L & svcf_status != "ok", NA_real_, final_svcf),

      expmt  = exper,
      sample = samp
    ) %>%
    ## Zero-reference rows that BAM evidence says are real. On a copy-neutral hemizygous locus
    ## SVCF = VAF, which at sv_ref = 0 is exactly 1: the tumour cells have lost their only copy, so
    ## the residual reads come from the normal fraction alone and observing no reference read is the
    ## physically correct outcome, not a technical failure. Recovering these is a scientific
    ## decision backed by the matched normal, which is why it arrives as a generated file rather
    ## than a code branch.
    { dat <- .
      if (!is.null(zero_ref_allowlist) && nrow(zero_ref_allowlist)) {
        al <- zero_ref_allowlist
        need <- c("sample", "chrom", "pos", "verdict")
        if (!all(need %in% names(al)))
          stop("calc_svcf: zero_ref_allowlist needs columns ",
               paste(need, collapse = ", "), call. = FALSE)
        al <- al[al$verdict == "recover", , drop = FALSE]
        if (nrow(al)) {
          key <- paste(dat$sample, dat$CHROM, dat$POS)
          hit <- key %in% paste(al$sample, al$chrom, al$pos)
          idx <- which(hit & dat$pl == 1L & dat$svcf_status == "zero_ref_depth")
          if (length(idx)) {
            dat$final_svcf[idx]  <- round(dat$sv_alt[idx] / (dat$sv_alt[idx] + dat$sv_ref[idx]), 2)
            dat$svcf_status[idx] <- "ok_zero_ref_recovered"
            message(sprintf("calc_svcf [%s]: recovered %d zero-reference chrX row(s) from the allowlist",
                            samp, length(idx)))
          }
        }
      }
      dat } %>%
    ## ------------------------------------------------------------------------------------------
    ## Hemizygous rows on a copy-altered segment. Three distinct cases, and conflating them is what
    ## went wrong before: a deletion or a duplication IS its own copy-number change, whereas an
    ## inversion or translocation merely sits inside one somebody else made.
    ## ------------------------------------------------------------------------------------------
    { dat <- .
      dat$sv_cnv_order  <- NA_character_
      dat$svcf_is_bound <- FALSE
      ## Hemizygous duplications, wherever their segment was classified: H3 needs cn_bar.
      i_hdup <- which(dat$pl == 1L & dat$classification == "DUP")
      if (length(i_hdup) && !is.null(hemi_cn_bar)) {
        cnb_d <- if (is.data.frame(hemi_cn_bar)) {
          hemi_cn_bar$cn_bar[match(paste(dat$CHROM[i_hdup], dat$POS[i_hdup]),
                                   paste(hemi_cn_bar$CHROM, hemi_cn_bar$POS))]
        } else rep_len(hemi_cn_bar, nrow(dat))[i_hdup]
        dd <- hemizygous_dup_svcf(cnb_d, hemi_dup_r)
        keep <- !is.na(cnb_d)
        dat$final_svcf[i_hdup][keep]    <- round(dd$svcf[keep], 2)
        dat$svcf_status[i_hdup][keep]   <- dd$status[keep]
        dat$sv_cnv_order[i_hdup][keep]  <- "dup_is_the_cnv"
        dat$svcf_is_bound[i_hdup][keep] <- dd$is_upper_bound[keep]
      }

      ## GATE ON MEASURED DEPTH, NOT ON cn_type. UPDATED 2026-07-31.
      ##
      ## This selected `cn_type != "norm"`, which is unreachable on a hemizygous chromosome. The
      ## SNP-based caller derives cn_type from heterozygous germline SNPs, and a male X has none,
      ## so it returns "norm" for every chrX row -- the ABSENCE of a call, not evidence of
      ## copy-neutrality. That is the same fact that makes ACR inestimable and is the reason this
      ## whole hemizygous path exists, so keying the gate on it made Eqs. 4a and 5a dead code
      ## exactly where they were needed.
      ##
      ## Measured on the COMBAT cohort: all 132 chrX rows are cn_type == "norm", while the depth
      ## says 111 of them have cn_bar > 1.05 (median 1.541, range 0.021-13.523) and only 16 sit
      ## within 0.05 of 1. Under the old gate all 119 non-DUP rows took SVCF = VAF.
      ##
      ## Safe as a superset: resolve_hemizygous_svcf() and hemizygous_del_svcf() each test
      ## copy-neutrality themselves and return VAF when |cn_bar - 1| <= tol, so a locus that
      ## really is copy-neutral gets the same answer it did before. Autosomes are untouched --
      ## the pl == 1L term is unchanged.
      idx <- which(dat$pl == 1L & dat$classification != "DUP")
      ## The rows the SNP-based caller DID flag. Only meaningful for the no-depth warning below,
      ## which used to describe exactly this set.
      idx_called <- which(dat$pl == 1L & !is.na(dat$cn_type) & dat$cn_type != "norm" &
                          dat$classification != "DUP")
      if (length(idx)) {
        if (is.null(hemi_cn_bar)) {
          if (length(idx_called)) {
            warning(sprintf(paste0("calc_svcf [%s]: %d hemizygous SV(s) lie on a copy-altered ",
                                   "segment but no hemi_cn_bar was supplied. They are left ",
                                   "unresolved rather than guessed. Run the chrX depth segmentation ",
                                   "and join cn_bar per SV; see JOIN-REVIEW.md."),
                            samp, length(idx_called)), call. = FALSE)
          }
        } else {
          ## hemi_cn_bar may be a per-row numeric vector, or the output of the segment-to-SV join:
          ## a table keyed on (CHROM, POS). The join is where "which segment applies when an SV
          ## spans a boundary" is decided (JOIN-REVIEW.md); this is only the lookup. An unmatched
          ## row is left unresolved and counted, never defaulted.
          if (is.data.frame(hemi_cn_bar)) {
            need <- c("CHROM", "POS", "cn_bar")
            if (!all(need %in% names(hemi_cn_bar))) {
              stop("calc_svcf: hemi_cn_bar table must have columns CHROM, POS, cn_bar; got: ",
                   paste(names(hemi_cn_bar), collapse = ", "), call. = FALSE)
            }
            key <- paste(dat$CHROM[idx], dat$POS[idx])
            cnb <- hemi_cn_bar$cn_bar[match(key, paste(hemi_cn_bar$CHROM, hemi_cn_bar$POS))]
            ## Warn only about rows the caller flagged as copy-altered. Since the gate widened to
            ## every hemizygous row, a missing cn_bar is now the ordinary case for a locus with no
            ## depth rather than a failure, and warning on all of them would bury the real signal.
            n_miss <- sum(is.na(cnb[match(idx_called, idx)]))
            if (n_miss) {
              warning(sprintf(paste0("calc_svcf [%s]: %d of %d hemizygous copy-altered SV(s) have ",
                                     "no cn_bar in the join table and are left unresolved."),
                              samp, n_miss, length(idx_called)), call. = FALSE)
            }
          } else {
            cnb <- rep_len(hemi_cn_bar, nrow(dat))[idx]
          }

          ## kappa for the CNV-first deletion form: the copies the locus would have WITHOUT the
          ## deletion, i.e. the flanking copy number. Resolved exactly as hemi_cn_bar is, and
          ## OPTIONAL -- absent, hemizygous_del_svcf() falls back to its previous behaviour, so a
          ## caller with no background estimate is unaffected.
          ##
          ## It matters wherever deletions land on copy-altered loci: on this cohort 13 of the 14
          ## chrX deletions do, and without kappa each of them falls through to h2, the form that
          ## scores 19.9% within 0.05 in the simulation against 58.9% for the CNV-first form.
          ## Produced by 08_chrx_segment_sv_join.py --bg-out.
          bgc <- NA_real_
          if (!is.null(hemi_bg_cn)) {
            if (is.data.frame(hemi_bg_cn)) {
              need_bg <- c("CHROM", "POS", "bg_cn")
              if (!all(need_bg %in% names(hemi_bg_cn))) {
                stop("calc_svcf: hemi_bg_cn table must have columns CHROM, POS, bg_cn; got: ",
                     paste(names(hemi_bg_cn), collapse = ", "), call. = FALSE)
              }
              key_bg <- paste(dat$CHROM[idx], dat$POS[idx])
              bgc <- hemi_bg_cn$bg_cn[match(key_bg, paste(hemi_bg_cn$CHROM, hemi_bg_cn$POS))]
            } else {
              bgc <- rep_len(hemi_bg_cn, nrow(dat))[idx]
            }
          }
          bgc <- rep_len(bgc, length(idx))

          ## Keep only rows with a measured cn_bar. A hemizygous row without one retains the
          ## copy-neutral SVCF = VAF assigned upstream rather than falling to "no_depth". Under
          ## the old gate such a row was never selected at all, so dropping it here would be a
          ## regression introduced by widening the gate rather than a decision about the data.
          ## With every branch below driven off `idx`, an empty idx is a no-op.
          keep_cnb <- is.finite(cnb)
          if (!all(keep_cnb)) {
            idx <- idx[keep_cnb]; cnb <- cnb[keep_cnb]; bgc <- bgc[keep_cnb]
          }

          cls <- dat$classification[idx]
          bpc <- dat$sv_alt[idx]; bec <- dat$sv_ref[idx]

          ## a hemizygous tandem duplication: SVCF = (cn_bar - 1)/(r - 1). r is not identifiable
          ## from one locus and r = 2 maximises SVCF, so these are upper bounds.
          i_dup <- which(cls == "DUP")
          if (length(i_dup)) {
            d <- hemizygous_dup_svcf(cnb[i_dup], hemi_dup_r)
            dat$final_svcf[idx[i_dup]]    <- round(d$svcf, 2)
            dat$svcf_status[idx[i_dup]]   <- d$status
            dat$sv_cnv_order[idx[i_dup]]  <- "dup_is_the_cnv"
            dat$svcf_is_bound[idx[i_dup]] <- d$is_upper_bound
          }

          ## a hemizygous deletion removes the only copy: reads and depth give SVCF independently,
          ## so their agreement is a genuine cross-check.
          i_del <- which(cls == "DEL")
          if (length(i_del)) {
            d <- hemizygous_del_svcf(bpc[i_del], bec[i_del], cnb[i_del], bg_cn = bgc[i_del])
            dat$final_svcf[idx[i_del]]   <- round(d$svcf, 2)
            dat$svcf_status[idx[i_del]]  <- d$status
            dat$sv_cnv_order[idx[i_del]] <- "del_is_the_cnv"
          }

          ## everything else sits inside a copy-number change made by another event. The sign of the
          ## SV-first form decides the ordering; no f_CNV and no integer copy number is needed.
          i_oth <- setdiff(seq_along(idx), c(i_dup, i_del))
          if (length(i_oth)) {
            res <- resolve_hemizygous_svcf(bpc[i_oth], bec[i_oth], cnb[i_oth])
            dat$final_svcf[idx[i_oth]]   <- round(res$svcf, 2)
            dat$svcf_status[idx[i_oth]]  <- res$status
            dat$sv_cnv_order[idx[i_oth]] <- res$ordering
          }
        }
      }
      dat } %>%
    ## infinities can only arise from ref == 0, which svcf_status already flags.
    ##
    ## Note this relabels on is.infinite(r_bar), so it fires even where final_svcf is finite and
    ## usable. On a copy-neutral hemizygous locus resolve_hemizygous_svcf() returns SVCF = VAF = 1
    ## with status "ok" at ref == 0, because VAF, unlike r_bar, is defined there. That result is
    ## deliberately suppressed. The reason is statistical rather than algebraic; see the MODELLING
    ## NOTE on svcf_status() in hemizygous.R before treating these rows as unrecoverable.
    mutate(
      ## Same containment: label everywhere, alter only hemizygous rows.
      ## This relabel fires on is.infinite(r_bar), which is true for EVERY sv_ref == 0 row,
      ## including ones the allowlist has just recovered whose final_svcf is finite and correct.
      ## Without the guard it silently undoes the recovery: the value survives but the status reads
      ## zero_ref_depth, so every downstream count treats the row as suppressed. Caught in test.
      svcf_status = ifelse(svcf_status != "ok_zero_ref_recovered" &
                             (is.infinite(final_svcf) | is.infinite(r_bar)),
                           "zero_ref_depth", svcf_status),
      final_svcf  = ifelse(pl == 1L & is.infinite(final_svcf), NA_real_, final_svcf)
    ) %>%
    ## SCOPE CONTAINMENT, second instance, and narrower than it first appears.
    ##
    ## The submitted beds RETAIN rows with r_bar = Inf (21 of 452 in 87955), so the original
    ## filter(!is.infinite(final_svcf), !is.infinite(r_bar)) was not in force when they were made.
    ## Reinstating it drops 203 autosomal rows and fails the gate on every sample. Do not.
    ##
    ## What the submitted beds do NOT retain is the one autosomal row whose raw_svcf is NaN: a
    ## tandem duplication with sv_ref = 0, where r_bar and r_2 are both Inf and raw_svcf = Inf/Inf.
    ## Exactly one such row exists in the cohort, chr4:10,075,098 in 87955, and the gate caught it.
    ## Drop it on autosomes to preserve the row count; keep it on hemizygous chromosomes, where a
    ## silently dropped row is what hid the chrX problem in the first place.
    filter(pl == 1L | !is.nan(raw_svcf)) %>%
    group_by(mate) %>%
    mutate(
      classification = ifelse(grepl('BND', classification), 'BND', classification),
      ## The original is mean(final_svcf) with NO na.rm, so a mate group containing an NA collapses
      ## to NA for the whole group. That behaviour must be preserved on autosomes. On hemizygous
      ## chromosomes we deliberately set NA for unresolvable rows (zero reference depth, duplications
      ## awaiting cn_bar), and letting one of those wipe out its partner's usable estimate would be
      ## a fresh error, so there na.rm = TRUE.
      final_svcf = ifelse(classification %in% c('INV', 'BND'),
                          ifelse(pl == 1L, mean(final_svcf, na.rm = TRUE), mean(final_svcf)),
                          final_svcf)
    ) %>%
    ungroup() %>%
    mutate(final_svcf = ifelse(is.nan(final_svcf), NA_real_, final_svcf))

  n_drop <- sum(final$svcf_status != "ok", na.rm = TRUE)
  if (n_drop > 0) {
    message(sprintf(
      "calc_svcf [%s]: %d of %d SV rows have no usable SVCF (%s)",
      samp, n_drop, nrow(final),
      paste(sprintf("%s=%d", names(table(final$svcf_status[final$svcf_status != "ok"])),
                    as.integer(table(final$svcf_status[final$svcf_status != "ok"]))),
            collapse = ", ")))
  }
  n_bound <- sum(final$svcf_is_bound, na.rm = TRUE)
  if (n_bound > 0) {
    message(sprintf("calc_svcf [%s]: %d hemizygous duplication(s) reported as UPPER BOUNDS (r=%g)",
                    samp, n_bound, hemi_dup_r))
  }

  return(final)
}
