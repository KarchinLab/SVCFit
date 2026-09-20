hemizygous_counts <- function(s, f, copies, order) {
  if (order == "sv_first") {
    stopifnot(f <= s)
    bpc <- copies * f + s - f
    bec <- 1 - s
  } else {
    stopifnot(s <= f)
    bpc <- s
    bec <- (copies - 1) * s + copies * (f - s) + 1 - f
  }
  list(bpc = bpc, bec = bec, cn_bar = 1 + f * (copies - 1))
}

test_that("hemizygous H1 and H2 recover counted cell mixtures", {
  cases <- data.frame(
    s = c(0.6, 0.4, 0.8, 0.2),
    f = c(0.4, 0.6, 0.5, 0.75),
    copies = c(3, 4, 2, 6),
    order = c("sv_first", "cnv_first", "sv_first", "cnv_first")
  )

  for (i in seq_len(nrow(cases))) {
    z <- cases[i, ]
    x <- hemizygous_counts(z$s, z$f, z$copies, z$order)
    out <- resolve_hemizygous_svcf(x$bpc, x$bec, x$cn_bar)
    expect_equal(x$bpc + x$bec, x$cn_bar)
    expect_equal(out$svcf, z$s)
    expect_equal(out$ordering,
                 if (z$order == "sv_first") "sv_before_cnv" else "cnv_before_sv")
  }
})

test_that("copy-neutral hemizygous SVCF reduces to VAF", {
  for (s in c(0.15, 0.5, 0.9, 1)) {
    out <- resolve_hemizygous_svcf(s, 1 - s, 1)
    expect_equal(out$svcf, s)
    expect_equal(out$ordering, "copy_neutral")
  }
})

test_that("hemizygous duplication uses carrier copy number and marks the bound", {
  for (copies in c(2, 3, 4, 6)) {
    for (s in c(0.1, 0.5, 1)) {
      cn_bar <- 1 + s * (copies - 1)
      exact <- hemizygous_dup_svcf(cn_bar, copies)
      upper <- hemizygous_dup_svcf(cn_bar, 2)
      expect_equal(exact$svcf, s)
      expect_true(upper$is_upper_bound)
      if (cn_bar <= 2) {
        expect_true(upper$svcf >= s - 1e-12)
        expect_equal(upper$status, "ok")
      } else {
        expect_true(is.na(upper$svcf))
        expect_equal(upper$status, "infeasible_svcf")
      }
    }
  }
})

test_that("hemizygous deletion uses depth and flags disagreement", {
  for (s in c(0.1, 0.4, 0.75, 0.95)) {
    out <- hemizygous_del_svcf(s, 1 - s, 1 - s)
    expect_equal(out$svcf, s)
    expect_equal(out$status, "ok")
  }

  mismatch <- hemizygous_del_svcf(0.5, 0.5, 0.9)
  expect_equal(mismatch$status, "del_depth_mismatch")
})
