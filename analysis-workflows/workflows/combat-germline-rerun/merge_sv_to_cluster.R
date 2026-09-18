#!/usr/bin/env Rscript
###############################################################################
# merge_sv2cluster.R
#
# Walk the SVCFit output base directory, find every pair's
# clustering/sv2cluster.csv, and concatenate them into a single CSV with an
# added pair_id column (set to the pair's directory name).
#
# Usage:
#   Rscript merge_sv2cluster.R \
#     --out_base /path/to/svcfit/output \
#     [--out_file <out_base>/sv2cluster_merged.csv]
###############################################################################
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(purrr)
})

option_list <- list(
  make_option(c("--out_base"), type = "character",
              default = Sys.getenv("GERMLINE_OUTPUT_DIR", unset = ""),
              help = "Base output dir containing per-pair subdirs; defaults to GERMLINE_OUTPUT_DIR"),
  make_option(c("--out_file"), type = "character", default = NULL,
              help = "Path to merged CSV [default: <out_base>/sv2cluster_merged.csv]")
)
opts <- parse_args(OptionParser(option_list = option_list))

if (!nzchar(opts$out_base)) stop("--out_base or GERMLINE_OUTPUT_DIR is required")

if (is.null(opts$out_file)) {
  opts$out_file <- file.path(opts$out_base, "sv2cluster_merged.csv")
}

# Find every clustering/sv2cluster.csv anywhere under out_base
csv_files <- list.files(
  path       = opts$out_base,
  pattern    = "^sv2cluster\\.csv$",
  recursive  = TRUE,
  full.names = TRUE
)
# Keep only files that live in a clustering/ folder (not stray copies)
csv_files <- csv_files[basename(dirname(csv_files)) == "clustering"]

if (length(csv_files) == 0L) {
  stop(sprintf("No clustering/sv2cluster.csv files found under %s", opts$out_base))
}

message(sprintf("Found %d sv2cluster.csv files to merge", length(csv_files)))

# Read each file, tag with pair_id (= the pair directory name)
# Force sample_ID to character — some pairs have purely numeric sample names
# (e.g. "81299"), which read_csv would otherwise infer as numeric and break
# bind_rows when other pairs have suffixed names like "81875_recut".
merged <- map_dfr(csv_files, function(f) {
  # Path looks like: <out_base>/<pair_id>/clustering/sv2cluster.csv
  pair_id <- basename(dirname(dirname(f)))
  df <- read_csv(f, show_col_types = FALSE,
                 col_types = cols(sample_ID = col_character(),
                                  .default  = col_guess()))
  message(sprintf("  %s: %d rows", pair_id, nrow(df)))
  df %>% mutate(pair_id = pair_id, .before = 1)
})

write_csv(merged, opts$out_file)
message(sprintf("\nMerged %d rows from %d pairs -> %s",
                nrow(merged), length(csv_files), opts$out_file))
