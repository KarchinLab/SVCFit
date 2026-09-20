test_that("resolve_svcfit_bed accepts a direct BED directory", {
  root <- tempfile("svcfit-direct-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  path <- file.path(root, "sample_a.bed")
  writeLines("CHROM\tPOS", path)

  expect_equal(SVCFit:::resolve_svcfit_bed(root, "sample_a"), normalizePath(path))
})

test_that("resolve_svcfit_bed accepts supported run-root layouts", {
  root <- tempfile("svcfit-nested-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nested <- file.path(root, "COMBAT", "SVCFit_output")
  dir.create(nested, recursive = TRUE)
  path <- file.path(nested, "sample_b.bed")
  writeLines("CHROM\tPOS", path)

  expect_equal(SVCFit:::resolve_svcfit_bed(root, "sample_b"), normalizePath(path))
})

test_that("resolve_svcfit_bed reports a missing sample clearly", {
  root <- tempfile("svcfit-missing-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  expect_error(SVCFit:::resolve_svcfit_bed(root, "absent"), "No SVCFit BED found")
})
