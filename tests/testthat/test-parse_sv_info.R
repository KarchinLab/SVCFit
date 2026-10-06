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

# One simulated VISOR translocation (MantaBND:31, rep1/exp1/c50p80m10): two junction
# records, each written on both mates, with RO/AO/RS/RP = 89/19/44/44 and 86/26/42/43.
read_bnd_fixture <- function() {
  sv <- read.table(test_path("fixtures", "svtyper_bnd.vcf"), quote = "\"")
  colnames(sv) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','normal','tumor')
  sv
}
parse_bnd <- function(sv) {
  bnd_info <- proc_bnd(sv, flank_del = 50)
  suppressMessages(parse_sv_info(sv, bnd_info$bnd, bnd_info$del, QUAL_thresh = 100, min_alt = 2))
}

test_that("parse_sv_info uses RS/2 + RP as the BND reference count", {
  sv_info <- parse_bnd(read_bnd_fixture())

  expect_equal(nrow(sv_info), 1)
  # selection keeps the higher-ref record (44/2 + 44 = 66, alt 19) and swaps in the
  # counts of the higher-alt record (42/2 + 43 = 64, alt 26)
  expect_equal(sv_info$sv_ref, 42 / 2 + 43)
  expect_equal(sv_info$sv_alt, 26)
})

test_that("parse_sv_info keeps RO for BND records without RS/RP", {
  sv <- read_bnd_fixture()
  drop <- function(fmt, val) {
    k <- strsplit(fmt, ":")[[1]]; v <- strsplit(val, ":")[[1]]
    paste(v[!k %in% c("RS", "RP")], collapse = ":")
  }
  sv$tumor  <- mapply(drop, sv$FORMAT, sv$tumor, USE.NAMES = FALSE)
  sv$normal <- mapply(drop, sv$FORMAT, sv$normal, USE.NAMES = FALSE)
  sv$FORMAT <- vapply(strsplit(sv$FORMAT, ":"), function(k) paste(k[!k %in% c("RS", "RP")], collapse = ":"), "")
  sv_info <- parse_bnd(sv)

  expect_equal(sv_info$sv_ref, 86)
  expect_equal(sv_info$sv_alt, 26)
})
