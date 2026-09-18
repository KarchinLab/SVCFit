library(dplyr)
library(tidyr)
library(stringr)
library(optparse)

option_list <- list(
  make_option(c("-s", "--surv"),
              type = "character",
              help = "Path to survivor output VCF"),

  make_option(c("-d", "--delly"),
              type = "character",
              help = "Path to delly output VCF"),

  make_option(c("-m", "--manta"),
              type = "character",
              help = "Path to manta output VCF"),

  make_option(c("-g", "--gridss"),
              type = "character",
              help = "Path to gridss output VCF"),

  make_option(c("-o", "--output"),
              type = "character",
              help = "Output VCF path"),

  make_option(c("-r", "--r_delly"),
              type = "character",
              help = "Path to delly VCF used for recovery"),

  make_option(c("-j", "--recovery_file"),
              type = "character",
              help = "Path to list of SVs to recover")
)

parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
opts <- parse_args(parser)

#-----------------------------
# helper functions
#-----------------------------
read_vcf_body <- function(vcf_file, col_names) {
  x <- read.table(
    vcf_file,
    sep = "\t",
    quote = "",
    comment.char = "#",
    header = FALSE,
    stringsAsFactors = FALSE,
    fill = TRUE
  )
  colnames(x) <- col_names
  x
}

get_vcf_header <- function(vcf_file, sample_name = "tumor") {
  hdr <- readLines(vcf_file)
  hdr <- hdr[grepl("^#", hdr)]

  meta_lines <- hdr[grepl("^##", hdr)]
  chrom_line <- hdr[grepl("^#CHROM", hdr)]

  # Replace the final column header line so it only has one sample column: tumor
  new_chrom_line <- paste(
    c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_name),
    collapse = "\t"
  )

  c(meta_lines, new_chrom_line)
}

write_vcf <- function(df, out_file, header_lines) {
  con <- file(out_file, open = "wt")
  on.exit(close(con))

  writeLines(header_lines, con)
  write.table(
    df,
    file = con,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

has_nonempty_file <- function(f) {
  if (is.null(f) || is.na(f) || !file.exists(f)) return(FALSE)
  length(readLines(f, warn = FALSE)) > 0
}

#-----------------------------
# read input VCFs
#-----------------------------
input <- read_vcf_body(
  opts$surv,
  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "delly", "gridss", "manta")
)

del <- read_vcf_body(
  opts$delly,
  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "tumor")
)

gri <- read_vcf_body(
  opts$gridss,
  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "tumor")
)

man <- read_vcf_body(
  opts$manta,
  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "tumor")
)

#-----------------------------
# keep header from survivor input
#-----------------------------
vcf_header <- get_vcf_header(opts$surv, sample_name = "tumor")

#-----------------------------
# get information from survivor output
#-----------------------------
tmp <- input %>%
  mutate(
    chr2   = gsub(".*CHR2=(.*?);END.*", "\\1", INFO),
    pos2   = gsub(".*END=(.*?);CIPOS.*", "\\1", INFO),
    svtype = gsub(".*TYPE=(.*?);SVM.*", "\\1", INFO),
    del_id = gsub("^(?:[^:]+:){7}([^:]+).*", "\\1", delly),
    gri_id = gsub("^(?:[^:]+:){7}([^:]+).*", "\\1", gridss),
    man_id = gsub("^(?:[^:]+:){7}([^:]+).*", "\\1", manta),
    man_id = gsub("_", ":", man_id)
  )

tmp_d <- del %>%
  filter(ID %in% tmp$del_id,
         !grepl("INS", ID)) %>%
  mutate(
    dRO = gsub("^(?:[^:]+:){14}([^:]+).*", "\\1", tumor),
    dRO = as.integer(dRO),
    dAO = gsub("^(?:[^:]+:){15}([^:]+).*", "\\1", tumor),
    dAO = as.integer(dAO),
    del_id = ID
  ) %>%
  filter(dAO >= 2)

tmp_g <- gri %>%
  filter(ID %in% tmp$gri_id,
         !grepl("INS", ID)) %>%
  mutate(
    gRO = gsub("^(?:[^:]+:){34}([^:]+).*", "\\1", tumor),
    gRO = as.integer(gRO),
    gAO = gsub("^(?:[^:]+:){35}([^:]+).*", "\\1", tumor),
    gAO = as.integer(gAO),
    gri_id = ID
  ) %>%
  filter(gAO >= 2)

tmp_m <- man %>%
  filter(ID %in% tmp$man_id,
         !grepl("INS", ID)) %>%
  mutate(
    mRO = ifelse(
      grepl("SR", FORMAT),
      gsub("^(?:[^:]+:){7}([^:]+).*", "\\1", tumor),
      gsub("^(?:[^:]+:){6}([^:]+).*", "\\1", tumor)
    ),
    mRO = as.integer(mRO),
    mAO = ifelse(
      grepl("SR", FORMAT),
      gsub("^(?:[^:]+:){8}([^:]+).*", "\\1", tumor),
      gsub("^(?:[^:]+:){7}([^:]+).*", "\\1", tumor)
    ),
    mAO = as.integer(mAO),
    man_id = ID
  ) %>%
  filter(mAO >= 2)

tmp2 <- tmp %>%
  left_join(tmp_d[, c("del_id", "dRO", "dAO")], by = "del_id") %>%
  left_join(tmp_g[, c("gri_id", "gRO", "gAO")], by = "gri_id") %>%
  left_join(tmp_m[, c("man_id", "mRO", "mAO")], by = "man_id") %>%
  mutate(
    RO = rowMeans(across(contains("RO"), ~ as.numeric(.x)), na.rm = TRUE),
    AO = rowMeans(across(contains("AO"), ~ as.numeric(.x)), na.rm = TRUE),
    RO = round(RO),
    AO = round(AO)
  ) %>%
  filter(!is.nan(AO), AO > 2)

#-----------------------------
# build output
#-----------------------------
output <- tmp2 %>%
  mutate(
    FORMAT = "GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB",
    tmp_str = "0/1:35:35.88:-6,-3,-12:19:15:3:15:3:0:0:0:15:3:0.17"
  ) %>%
  rowwise() %>%
  mutate(
    tumor = {
      parts <- str_split(tmp_str, ":", simplify = TRUE)
      parts[5] <- as.character(RO + AO)  # DP
      parts[6] <- as.character(RO)       # RO
      parts[7] <- as.character(AO)       # AO
      str_c(parts, collapse = ":")
    }
  ) %>%
  ungroup() %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, tumor)

#-----------------------------
# recover SVs if requested
#-----------------------------
if (has_nonempty_file(opts$recovery_file)) {
  r_del <- read_vcf_body(
    opts$r_delly,
    c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "tumor", "normal")
  )

  recovery <- read.delim(opts$recovery_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(recovery) <- c("CHROM", "POS", "END", "SVtype")

  recovery <- recovery %>%
    mutate(key = as.character(floor(POS * 0.01)))

  recr_sv <- r_del %>%
    filter(
      CHROM %in% recovery$CHROM,
      gsub(".*TYPE=(.*);SVM.*", "\\1", INFO) %in% toupper(recovery$SVtype),
      as.character(floor(POS * 0.01)) %in% recovery$key
    ) %>%
    separate(
      tumor,
      into = c("GT", "GL", "GQ", "FT", "RCL", "RC", "RCR", "RDCN", "DR", "DV", "RR", "RV"),
      sep = ":",
      remove = FALSE
    ) %>%
    mutate(
      RO = as.numeric(DR) + as.numeric(RR),
      AO = as.numeric(DV) + as.numeric(RV)
    )

  add <- recr_sv %>%
    mutate(
      FORMAT = "GT:GQ:SQ:GL:DP:RO:AO:QR:QA:RS:AS:ASC:RP:AP:AB",
      tmp_str = "0/1:35:35.88:-6,-3,-12:19:15:3:15:3:0:0:0:15:3:0.17"
    ) %>%
    rowwise() %>%
    mutate(
      tumor = {
        parts <- str_split(tmp_str, ":", simplify = TRUE)
        parts[5] <- as.character(RO + AO)
        parts[6] <- as.character(RO)
        parts[7] <- as.character(AO)
        str_c(parts, collapse = ":")
      }
    ) %>%
    ungroup() %>%
    select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, tumor) %>%
    mutate(FILTER = "PASS")

  new_out <- bind_rows(output, add)
} else {
  message("no sv to recover.")
  new_out <- output
}

#-----------------------------
# write output VCF with header
#-----------------------------
write_vcf(new_out, opts$output, vcf_header)
