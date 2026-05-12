p_het  <- system.file("extdata", "examples/het_near_sv_c50p80m50.vcf", package = "SVCFit")
p_onsv <- system.file("extdata", "examples/het_on_sv_c50p80m50.vcf",   package = "SVCFit")
p_sv   <- system.file("extdata", "examples/example_sv.bed",             package = "SVCFit")
p_cnv  <- system.file("extdata", "examples/c50p80m50.bed",              package = "SVCFit")

test_that("proc_bnd returns a named list with bnd and del elements", {
  skip_if(nchar(p_sv) == 0, "example data not installed")

  data   <- load_data(p_het, p_onsv, p_sv, p_cnv, chr = NULL, tumor_only = FALSE)
  result <- proc_bnd(data$sv, flank_del = 50)

  expect_type(result, "list")
  expect_named(result, c("bnd", "del"))
  expect_s3_class(result$bnd, "data.frame")
  expect_s3_class(result$del, "data.frame")
})

test_that("proc_bnd bnd records contain required columns", {
  skip_if(nchar(p_sv) == 0, "example data not installed")

  data   <- load_data(p_het, p_onsv, p_sv, p_cnv, chr = NULL, tumor_only = FALSE)
  result <- proc_bnd(data$sv, flank_del = 50)

  if (nrow(result$bnd) > 0) {
    expect_true(all(c("CHROM", "POS", "ID", "class") %in% names(result$bnd)))
    expect_true(all(result$bnd$class %in% c("rtl", "coptl", "cuttl")))
  } else {
    skip("no BND records in example data")
  }
})
