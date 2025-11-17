#' Title
#'
#' @param sv an object of class 'data frame' which stores the sv input from `load_data`.
#' @param bnd_thresh an object of class 'integer' which describes.
#'
#' @returns sth
#' @export
#'
preproc_bnd <- function(sv) {
  
  tmp=sv %>%
    filter(grepl("BND", INFO)) %>%
    mutate(chr2=paste0('chr',gsub('.*chr(.*):.*','\\1', ALT)),
           pos2=gsub('.*chr.*:(\\d+).*','\\1', ALT),
           pos2=as.integer(pos2))
  #chr2=gsub('.*CHR2=(.*);END.*','\\1', INFO),
  #pos2=gsub('.*END=(.*);CIPOS.*','\\1', INFO),
  #pos2=as.integer(pos2),
  
  bnd_tmp=tmp %>%
    rowwise()%>%
    mutate(rn1=which(CHROM==tmp$CHROM & abs(POS-tmp$POS)<400 & chr2==tmp$chr2)[1]*2-1,
           rn2=which(chr2==tmp$chr2 & abs(pos2-tmp$pos2)<400 & CHROM==tmp$CHROM)[1]*2-1)%>%
    group_by(rn2)%>%
    mutate(len1=max(POS)-min(POS))%>%
    group_by(rn1)%>%
    mutate(len2=max(pos2)-min(pos2),
           sum_rn2=sum(as.integer(rn2))) %>%
    ungroup()%>%
    mutate(mateid=gsub(".*MATEID=([^;]+);.*", "\\1", INFO))%>%
    rowwise()%>%
    mutate(t_sum_rn2=.$sum_rn2[which(mateid==.$ID)[1]],
           n_sum_rn2=sum(sum_rn2, t_sum_rn2))%>%
    group_by(n_sum_rn2)%>%
    mutate(donor=unique(CHROM[which.max(len1)]),
           receiver=unique(CHROM[which.min(len1)]),
           receiver=ifelse(receiver==donor, unique(chr2[chr2!=receiver]), receiver))%>%
    # donor=unique(CHROM[which.max(len1)]),
    # receiver=unique(CHROM[which.min(len1)]),
    # receiver=ifelse(class=='rtl', CHROM[which(CHROM!=donor)], receiver))%>%
    filter(CHROM==receiver)%>%
    mutate(n=n(),
           class=ifelse(n>2, 'rtl','ubtl'),
           class=ifelse(max(len2)<50, 'rtl',class))
  
  return(bnd_tmp)
}

# 2b. Identify DEL overlapping BND
#' Title
#'
#' @param sv an object of class 'data frame' which stores the sv input from `load_data`.
#' @param bnd_tmp an object of class 'data frame' which is the output from `preproc_bnd`.
#' @param flank_del an object of class 'integer' which describes the maximum genomic
#' location difference between a deletion and translocation to be considered related.
#'
#' @returns sth
#' @export
#'
process_del <- function(sv, bnd_tmp, flank = 50) {
  
  del=sv %>%
    filter(str_detect(INFO, 'DEL'))%>%
    mutate(chr2=CHROM,
           pos2=gsub('.*END=(.*);CIPOS.*','\\1',INFO),
           pos2=as.integer(pos2))%>%
    rowwise() %>%
    mutate(
      #check any del is associated with trans (cut-paste)
      overlap = any(
        CHROM == bnd_tmp$chr2 &
          abs(POS-bnd_tmp$pos2)<50
      )
    )%>%
    filter(overlap) %>%
    # assign trans id to associated del
    mutate(
      BNDid = bnd_tmp$ID[which(
        CHROM == bnd_tmp$chr2 &
          abs(POS-bnd_tmp$pos2)<50
      )]
    ) %>%
    ungroup()
  return(del)
}

# 2c. Annotate BND classes
#' Title
#'
#' @param bnd_tmp an object of class 'data frame' which stores the sv input from `preproc_bnd`.
#' @param del an object of class 'data frame' which stores the sv input from `process_del`.
#'
#' @returns sth
#' @export
#'
annotate_bnd <- function(bnd_tmp, del) {
  bnd=bnd_tmp %>%
    # use del to classify cut-paste, and rest to copy paste
    mutate(
      class = case_when(
        ID %in% del$BNDid ~ 'cuttl',
        class == 'ubtl'      ~ 'coptl',
        TRUE                 ~ class
      ))%>%
    group_by(n_sum_rn2)%>%
    mutate(class=ifelse(any(grepl('cuttl', class)), 'cuttl',class))%>%
    # figure out the segment length that was translocated
    # group_by(class) %>%
    # mutate(BND_num = as.integer(str_extract(ID, '\\d+')))%>%
    # arrange(BND_num, .by_group = TRUE) %>%
    # mutate(group = cumsum(c(0, diff(BND_num) > 1)) + 1) %>%
    # ungroup()
    group_by(class, n_sum_rn2, CHROM) %>%
    mutate(nPOS=min(pos2),
           mPOS=max(pos2))
  return(bnd)
}
