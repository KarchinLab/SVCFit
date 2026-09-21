#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || any(args %in% c("-h", "--help"))) {
  cat("Usage: compare_svcf_runs.R OLD_INPUT NEW_INPUT OUTPUT_DIR\n",
      "INPUT may be an RDS data.frame, a tab-delimited text file, or a run root\n",
      "containing per-sample BED files.\n", sep = "")
  quit(status = if (length(args) == 3L) 0L else 2L)
}

old_path <- normalizePath(args[[1]], mustWork = TRUE)
new_path <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Could not create output directory: ", out_dir, call. = FALSE)

read_one_table <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    readRDS(path)
  } else {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
}

bind_tables <- function(tables) {
  columns <- unique(unlist(lapply(tables, names), use.names = FALSE))
  tables <- lapply(tables, function(x) {
    missing <- setdiff(columns, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[columns]
  })
  do.call(rbind, tables)
}

read_input <- function(path, label) {
  if (dir.exists(path)) {
    candidates <- list.files(path, pattern = "\\.bed$", recursive = TRUE,
                             full.names = TRUE)
    normalized <- gsub("\\\\", "/", candidates)
    preferred <- candidates[grepl("/SVCFit_output/[^/]+\\.bed$", normalized)]
    files <- if (length(preferred)) preferred else candidates
    files <- sort(normalizePath(files, mustWork = TRUE))
    if (!length(files)) {
      stop("No per-sample BED files under ", path, call. = FALSE)
    }
    sample_names <- tools::file_path_sans_ext(basename(files))
    if (anyDuplicated(sample_names)) {
      stop(label, " run contains a sample more than once: ",
           paste(unique(sample_names[duplicated(sample_names)]), collapse = ", "), call. = FALSE)
    }
    tables <- lapply(seq_along(files), function(i) {
      x <- read_one_table(files[[i]])
      if (!"sample" %in% names(x)) x$sample <- sample_names[[i]]
      if (any(as.character(x$sample) != sample_names[[i]])) {
        stop("Sample column disagrees with file name: ", files[[i]], call. = FALSE)
      }
      x[[paste0("source_file_", label)]] <- files[[i]]
      x
    })
    x <- bind_tables(tables)
  } else {
    files <- path
    x <- read_one_table(path)
  }
  if (!is.data.frame(x)) stop("Input is not a data.frame: ", path, call. = FALSE)
  if (!"final_svcf" %in% names(x)) stop("Input lacks final_svcf: ", path, call. = FALSE)
  list(data = x, files = files)
}

old_input <- read_input(old_path, "old")
new_input <- read_input(new_path, "new")
old <- old_input$data
new <- new_input$data

write_manifest <- function(files, label) {
  info <- file.info(files)
  manifest <- data.frame(
    run = label,
    file = files,
    bytes = info$size,
    md5 = unname(tools::md5sum(files)),
    stringsAsFactors = FALSE
  )
  write.table(manifest, file.path(out_dir, paste0(label, "_input_manifest.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_manifest(old_input$files, "old")
write_manifest(new_input$files, "new")

write.table(old, file.path(out_dir, "old_combined.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")
write.table(new, file.path(out_dir, "new_combined.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

validity_metrics <- function(x, arm) {
  value <- suppressWarnings(as.numeric(x$final_svcf))
  finite <- is.finite(value)
  data.frame(
    arm = arm,
    rows = length(value),
    finite_final_svcf = sum(finite),
    missing_or_nonfinite_final_svcf = sum(!finite),
    below_zero = sum(finite & value < 0),
    above_one = sum(finite & value > 1),
    zero_ref_depth = if ("svcf_status" %in% names(x))
      sum(x$svcf_status == "zero_ref_depth", na.rm = TRUE) else NA_integer_,
    stringsAsFactors = FALSE
  )
}

count_column <- function(x, arm, column) {
  if (!column %in% names(x)) return(NULL)
  value <- as.character(x[[column]])
  value[is.na(value) | value == ""] <- "NA"
  counts <- as.data.frame(table(value), stringsAsFactors = FALSE)
  names(counts) <- c("status", "rows")
  counts$arm <- arm
  counts$column <- column
  counts[c("arm", "column", "status", "rows")]
}

validity <- rbind(validity_metrics(old, "old"), validity_metrics(new, "new"))
write.table(validity, file.path(out_dir, "svcf_old_vs_new_validity.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

status_counts <- do.call(rbind, Filter(Negate(is.null), list(
  count_column(old, "old", "svcf_status"),
  count_column(new, "new", "svcf_status"),
  count_column(old, "old", "final_svcf_constraint_status"),
  count_column(new, "new", "final_svcf_constraint_status")
)))
write.table(status_counts, file.path(out_dir, "svcf_old_vs_new_status_counts.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

preferred_keys <- c("sample", "expmt", "CHROM", "POS", "END", "ID", "mate")
keys <- preferred_keys[preferred_keys %in% names(old) & preferred_keys %in% names(new)]
if (!length(keys)) stop("The inputs have no common event-key columns", call. = FALSE)
if (anyDuplicated(old[keys]) || anyDuplicated(new[keys])) {
  stop("Common event keys are not unique: ", paste(keys, collapse = ", "), call. = FALSE)
}

old$.present_old <- TRUE
new$.present_new <- TRUE
joined <- merge(old, new, by = keys, all = TRUE, suffixes = c("_old", "_new"), sort = FALSE)
joined$match_status <- ifelse(is.na(joined$.present_old), "new_only",
                              ifelse(is.na(joined$.present_new), "old_only", "matched"))
old_value <- joined$final_svcf_old
new_value <- joined$final_svcf_new
joined$svcf_delta <- new_value - old_value
joined$estimate_changed <- joined$match_status != "matched" |
  xor(is.na(old_value), is.na(new_value)) |
  (!is.na(old_value) & !is.na(new_value) & abs(joined$svcf_delta) > 1e-12)

write.table(joined, file.path(out_dir, "svcf_old_vs_new_events.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

matched <- joined$match_status == "matched"
finite_delta <- matched & is.finite(joined$svcf_delta)
constraint_column <- if ("ss2_constraint_status_new" %in% names(joined)) {
  "ss2_constraint_status_new"
} else if ("ss2_constraint_status" %in% names(joined) &&
           !"ss2_constraint_status" %in% names(old)) {
  "ss2_constraint_status"
} else {
  NA_character_
}
constraint_status <- if (!is.na(constraint_column)) joined[[constraint_column]] else NULL

summary <- data.frame(
  old_rows = nrow(old), new_rows = nrow(new), matched_rows = sum(matched),
  old_only = sum(joined$match_status == "old_only"),
  new_only = sum(joined$match_status == "new_only"),
  changed_rows = sum(joined$estimate_changed),
  old_evaluable = sum(is.finite(old$final_svcf)),
  new_evaluable = sum(is.finite(new$final_svcf)),
  mean_absolute_change = if (any(finite_delta)) mean(abs(joined$svcf_delta[finite_delta])) else NA_real_,
  maximum_absolute_change = if (any(finite_delta)) max(abs(joined$svcf_delta[finite_delta])) else NA_real_,
  new_boundary_low = if (!is.null(constraint_status))
    sum(constraint_status == "boundary_low", na.rm = TRUE) else NA_integer_,
  new_boundary_high = if (!is.null(constraint_status))
    sum(constraint_status == "boundary_high", na.rm = TRUE) else NA_integer_,
  old_below_zero = validity$below_zero[validity$arm == "old"],
  old_above_one = validity$above_one[validity$arm == "old"],
  new_below_zero = validity$below_zero[validity$arm == "new"],
  new_above_one = validity$above_one[validity$arm == "new"],
  old_zero_ref_depth = validity$zero_ref_depth[validity$arm == "old"],
  new_zero_ref_depth = validity$zero_ref_depth[validity$arm == "new"]
)
write.table(summary, file.path(out_dir, "svcf_old_vs_new_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")

group_candidates <- c("classification", "cn_type", "bkg_cnv", "zygosity")
group_columns <- character()
for (nm in group_candidates) {
  if (nm %in% keys) group_columns <- c(group_columns, nm)
  if (paste0(nm, "_new") %in% names(joined)) group_columns <- c(group_columns, paste0(nm, "_new"))
}
if (length(group_columns)) {
  group_key <- interaction(joined[group_columns], drop = TRUE, lex.order = TRUE)
  by_group <- do.call(rbind, lapply(split(seq_len(nrow(joined)), group_key), function(i) {
    z <- joined[i, , drop = FALSE]
    d <- z$svcf_delta[is.finite(z$svcf_delta)]
    out <- z[1, group_columns, drop = FALSE]
    out$rows <- nrow(z)
    out$changed_rows <- sum(z$estimate_changed)
    out$mean_absolute_change <- if (length(d)) mean(abs(d)) else NA_real_
    out$maximum_absolute_change <- if (length(d)) max(abs(d)) else NA_real_
    out
  }))
  rownames(by_group) <- NULL
  write.table(by_group, file.path(out_dir, "svcf_old_vs_new_by_group.tsv"), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
}

report <- c(
  "# SVCF old-versus-new run comparison", "",
  paste0("- Old input: `", old_path, "`"),
  paste0("- New input: `", new_path, "`"),
  paste0("- Old input files: ", length(old_input$files),
         " (checksums in `old_input_manifest.tsv`)"),
  paste0("- New input files: ", length(new_input$files),
         " (checksums in `new_input_manifest.tsv`)"), "",
  "| Measure | Value |", "|---|---:|",
  paste0("| Old rows | ", summary$old_rows, " |"),
  paste0("| New rows | ", summary$new_rows, " |"),
  paste0("| Matched rows | ", summary$matched_rows, " |"),
  paste0("| Changed rows | ", summary$changed_rows, " |"),
  paste0("| Old finite values below 0 | ", summary$old_below_zero, " |"),
  paste0("| Old finite values above 1 | ", summary$old_above_one, " |"),
  paste0("| New finite values below 0 | ", summary$new_below_zero, " |"),
  paste0("| New finite values above 1 | ", summary$new_above_one, " |"),
  paste0("| Old zero-reference rows | ", summary$old_zero_ref_depth, " |"),
  paste0("| New zero-reference rows | ", summary$new_zero_ref_depth, " |"),
  paste0("| New upper-boundary estimates | ", summary$new_boundary_high, " |"),
  paste0("| Mean absolute SVCF change | ", signif(summary$mean_absolute_change, 6), " |"),
  paste0("| Maximum absolute SVCF change | ", signif(summary$maximum_absolute_change, 6), " |"), "",
  "Review the event table before replacing any accepted downstream output."
)
writeLines(report, file.path(out_dir, "SVCF_SHADOW_COMPARISON.md"))

cat("Wrote SVCF comparison to ", out_dir, "\n", sep = "")

if (!"final_svcf_constraint_status" %in% names(new)) {
  stop("New run lacks final_svcf_constraint_status; it was not produced by the globally constrained SVCFit implementation",
       call. = FALSE)
}
new_invalid <- summary$new_below_zero + summary$new_above_one
if (new_invalid > 0) {
  stop("New run contains ", new_invalid,
       " finite final_svcf value(s) outside [0, 1]; refusing shadow acceptance",
       call. = FALSE)
}
