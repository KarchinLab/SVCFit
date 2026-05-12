p_het  <- system.file("extdata", "examples/het_near_sv_c50p80m50.vcf", package = "SVCFit")
p_onsv <- system.file("extdata", "examples/het_on_sv_c50p80m50.vcf",   package = "SVCFit")
p_sv   <- system.file("extdata", "examples/example_sv.bed",             package = "SVCFit")
p_cnv  <- system.file("extdata", "examples/c50p80m50.bed",              package = "SVCFit")

setup_sv_data <- function() {
  data     <- load_data(p_het, p_onsv, p_sv, p_cnv, chr = NULL, tumor_only = FALSE)
  bnd_info <- proc_bnd(data$sv, flank_del = 50)
  list(data = data, bnd_info = bnd_info)
}

test_that("parse_sv_info returns a data.frame with required columns", {
  skip_if(nchar(p_sv) == 0, "example data not installed")

  d       <- setup_sv_data()
  sv_info <- parse_sv_info(d$data$sv, d$bnd_info$bnd, d$bnd_info$del,
                           QUAL_thresh = 100, min_alt = 2)

  expect_s3_class(sv_info, "data.frame")
  expect_true(all(c("CHROM", "POS", "END", "ID", "sv_ref", "sv_alt",
                    "classification") %in% names(sv_info)))
})

test_that("parse_sv_info min_alt filter removes low-support SVs", {
  skip_if(nchar(p_sv) == 0, "example data not installed")

  d        <- setup_sv_data()
  sv_loose <- parse_sv_info(d$data$sv, d$bnd_info$bnd, d$bnd_info$del,
                            QUAL_thresh = 100, min_alt = 2)
  sv_strict <- parse_sv_info(d$data$sv, d$bnd_info$bnd, d$bnd_info$del,
                             QUAL_thresh = 100, min_alt = 50)

  expect_gte(nrow(sv_loose), nrow(sv_strict))
  expect_true(all(sv_strict$sv_alt >= 50))
})

test_that("parse_sv_info excludes insertion SVs", {
  skip_if(nchar(p_sv) == 0, "example data not installed")

  d       <- setup_sv_data()
  sv_info <- parse_sv_info(d$data$sv, d$bnd_info$bnd, d$bnd_info$del,
                           QUAL_thresh = 100, min_alt = 2)

  expect_false(any(grepl("INS", sv_info$classification)))
})
