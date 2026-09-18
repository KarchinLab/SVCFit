#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(readr)
})

options <- list(
  make_option(c("-s", "--sample"), type = "character", help = "Sample name"),
  make_option(c("-o", "--o_path"), type = "character", help = "Output directory"),
  make_option(c("-v", "--sv_file"), type = "character", help = "Input SV VCF")
)
opts <- parse_args(OptionParser(option_list = options))

required <- c("sample", "o_path", "sv_file")
missing <- required[vapply(required, function(x) is.null(opts[[x]]) || !nzchar(opts[[x]]), logical(1))]
if (length(missing)) stop("Missing required options: ", paste(missing, collapse = ", "))
if (!file.exists(opts$sv_file)) stop("SV VCF does not exist: ", opts$sv_file)
dir.create(opts$o_path, recursive = TRUE, showWarnings = FALSE)

sv <- read.delim(
  opts$sv_file, header = FALSE, comment.char = "#", quote = "",
  stringsAsFactors = FALSE, fill = TRUE, check.names = FALSE
)
if (ncol(sv) < 8L) stop("Expected at least eight VCF columns in ", opts$sv_file)
names(sv)[1:8] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

is_bnd <- grepl("BND", sv$ID) | grepl("[\\[\\]]", sv$ALT)
mate_match <- regexec("[\\[\\]]([^:\\[\\]]+):(\\d+)[\\[\\]]", sv$ALT)
mate_parts <- regmatches(sv$ALT, mate_match)
mate_chr <- vapply(mate_parts, function(x) if (length(x) >= 3L) x[2] else NA_character_, character(1))
mate_pos <- vapply(mate_parts, function(x) if (length(x) >= 3L) x[3] else NA_character_, character(1))
end_match <- regexec("(?:^|;)END=([^;]+)", sv$INFO, perl = TRUE)
end_parts <- regmatches(sv$INFO, end_match)
end_pos <- vapply(end_parts, function(x) if (length(x) >= 2L) x[2] else NA_character_, character(1))

right_chr <- ifelse(is_bnd, mate_chr, sv$CHROM)
right_pos <- suppressWarnings(as.integer(ifelse(is_bnd, mate_pos, end_pos)))
left_pos <- suppressWarnings(as.integer(sv$POS))
if (anyNA(left_pos) || anyNA(right_pos) || anyNA(right_chr)) {
  stop("Could not parse one or more SV breakpoint coordinates in ", opts$sv_file)
}

# Preserve the legacy window convention: coordinate minus/plus 500, with starts
# clamped to zero. bcftools interprets the three-column output as BED regions.
regions <- rbind(
  data.frame(CHROM = sv$CHROM, lpos = pmax(0L, left_pos - 500L), rpos = left_pos + 500L),
  data.frame(CHROM = right_chr, lpos = pmax(0L, right_pos - 500L), rpos = right_pos + 500L)
)
regions <- unique(regions)
output <- file.path(opts$o_path, paste0(opts$sample, ".bed"))
write_delim(regions, output, delim = "\t", col_names = FALSE)
message("Wrote ", nrow(regions), " SV-flanking regions to ", output)
