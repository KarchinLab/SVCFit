p_het  <- system.file("extdata", "examples/het_near_sv_c50p80m50.vcf", package = "SVCFit")
p_onsv <- system.file("extdata", "examples/het_on_sv_c50p80m50.vcf",   package = "SVCFit")
p_sv   <- system.file("extdata", "examples/example_sv.bed",             package = "SVCFit")
p_cnv  <- system.file("extdata", "examples/c50p80m50.bed",              package = "SVCFit")

test_that("load_data returns a named list with four elements", {
  skip_if(nchar(p_sv) == 0, "example data not installed")

  data <- load_data(p_het, p_onsv, p_sv, p_cnv, chr = NULL, tumor_only = FALSE)

  expect_type(data, "list")
  expect_named(data, c("het_snp", "het_on_sv", "sv", "cnv"))
  expect_s3_class(data$het_snp,  "data.frame")
  expect_s3_class(data$het_on_sv,"data.frame")
  expect_s3_class(data$sv,       "data.frame")
  expect_s3_class(data$cnv,      "data.frame")
})

test_that("load_data SV columns are complete and chr-prefixed", {
  skip_if(nchar(p_sv) == 0, "example data not installed")

  data <- load_data(p_het, p_onsv, p_sv, p_cnv, chr = NULL, tumor_only = FALSE)

  expect_true(all(c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO",
                    "FORMAT","normal","tumor") %in% names(data$sv)))
  expect_true(all(grepl("^chr", data$sv$CHROM)))
  expect_true(all(grepl("^chr", data$cnv$chrom)))
})

test_that("load_data chr filter restricts output to requested chromosome", {
  skip_if(nchar(p_sv) == 0, "example data not installed")

  data <- load_data(p_het, p_onsv, p_sv, p_cnv, chr = "chr1", tumor_only = FALSE)

  expect_true(all(data$sv$CHROM      == "chr1"))
  expect_true(all(data$het_snp$CHROM == "chr1"))
  expect_true(all(data$cnv$chrom     == "chr1"))
})

test_that("load_data returns non-empty data frames for the example dataset", {
  skip_if(nchar(p_sv) == 0, "example data not installed")

  data <- load_data(p_het, p_onsv, p_sv, p_cnv, chr = NULL, tumor_only = FALSE)

  expect_gt(nrow(data$sv),      0)
  expect_gt(nrow(data$cnv),     0)
  expect_gt(nrow(data$het_snp), 0)
})
