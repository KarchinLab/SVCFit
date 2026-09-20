test_that("SVCF constraint projects finite candidates onto the unit interval", {
  x <- c(-0.2, 0, 0.25, 1, 1.2, 2, Inf, NaN, NA_real_)

  expect_equal(
    SVCFit:::constrain_svcf(x),
    c(0, 0, 0.25, 1, 1, 1, NA, NA, NA)
  )
  expect_equal(
    SVCFit:::svcf_constraint_status(x),
    c("boundary_low", "in_range", "in_range", "in_range",
      "boundary_high", "boundary_high", "nonfinite", "nonfinite", "nonfinite")
  )
})

test_that("zygosity correction is applied before the SVCF constraint", {
  s2_raw <- 2

  expect_equal(SVCFit:::constrain_svcf(s2_raw), 1)
  expect_equal(SVCFit:::constrain_svcf(0.5 * s2_raw), 1)

  # The retired modulo-first calculation returned zero in both cases.
  expect_equal(s2_raw %% 2, 0)
  expect_equal(0.5 * (s2_raw %% 2), 0)
})

test_that("calc_svcf retains raw s2 and flags an active upper boundary", {
  anno <- data.frame(
    CHROM = c("chr1", "chr1"), POS = c(100L, 200L), ID = c("het", "hom"),
    zygosity = c("het", "hom"), ASCN = c(1.5, 1.5), cn_type = c("DEL", "DEL"),
    bkg_cnv = c("dup", "dup"), major = c(1, 1), minor = c(0, 0),
    mate = c("het", "hom"), stringsAsFactors = FALSE
  )
  reads <- data.frame(
    CHROM = c("chr1", "chr1"), POS = c(100L, 200L), ID = c("het", "hom"),
    sv_alt = c(8, 8), sv_ref = c(2, 2), classification = c("DEL", "DEL"),
    stringsAsFactors = FALSE
  )

  out <- calc_svcf(anno, reads, samp = "sample", exper = "experiment")

  expect_equal(out$s2_raw, c(2, 2))
  expect_equal(out$ss2_raw, c(2, 1))
  expect_equal(out$ss2, c(1, 1))
  expect_equal(out$ss2_constraint_status, c("boundary_high", "in_range"))
  expect_equal(out$final_svcf, c(1, 1))
})

test_that("calc_svcf uses the sign criterion with the 0.1 noise buffer", {
  anno <- data.frame(
    CHROM = c("chr1", "chr1", "chr1"), POS = c(100L, 200L, 300L),
    ID = c("buffer", "positive", "deletion"), zygosity = "het",
    ASCN = c(1.79, 1.5, 1.2), cn_type = c("DUP", "DUP", "DEL"),
    bkg_cnv = "dup", major = c(3, 3, 1), minor = 0,
    mate = c("buffer", "positive", "deletion"), stringsAsFactors = FALSE
  )
  reads <- data.frame(
    CHROM = "chr1", POS = c(100L, 200L, 300L),
    ID = c("buffer", "positive", "deletion"),
    sv_alt = c(3, 4, 4), sv_ref = c(7, 6, 6),
    classification = c("INV", "INV", "DEL"), stringsAsFactors = FALSE
  )

  out <- calc_svcf(anno, reads, samp = "sample", exper = "experiment")
  out <- out[match(c("buffer", "positive", "deletion"), out$ID), ]

  expect_equal(out$ss1, c(0.047, 0.5, 0.68), tolerance = 1e-12)
  expect_equal(out$ss2, c(0.837, 1, 0.88), tolerance = 1e-12)
  expect_equal(out$final_svcf, c(0.84, 0.5, 0.88))
})

test_that("diploid duplication uses FACETS carrier copy number", {
  anno <- data.frame(
    CHROM = "chr1", POS = c(100L, 200L, 300L),
    ID = c("r3_clonal", "r3_half", "r4_half"), zygosity = "het",
    ASCN = 1, cn_type = "norm", bkg_cnv = "norm",
    major = c(3, 3, 4), minor = 0,
    mate = c("r3_clonal", "r3_half", "r4_half"), stringsAsFactors = FALSE
  )
  reads <- data.frame(
    CHROM = "chr1", POS = c(100L, 200L, 300L),
    ID = c("r3_clonal", "r3_half", "r4_half"),
    sv_alt = c(2, 1, 2), sv_ref = c(4, 4, 4),
    classification = "DUP", stringsAsFactors = FALSE
  )

  out <- calc_svcf(anno, reads, samp = "sample", exper = "experiment")
  out <- out[match(c("r3_clonal", "r3_half", "r4_half"), out$ID), ]

  expect_equal(out$r_2, c(1, 1, 2))
  expect_equal(out$final_svcf, c(1, 0.5, 0.5))
})
