#' Process translocation
#'
#' @param sv an object of class 'data frame' which stores the sv input from `load_data`.
#' @param flank_del an object of class 'integer' which describes the maximum genomic location 
#' difference between a deletion and a translocation to be considered overlapping 
#'
#' @returns sth
#' @export
#'
proc_bnd <- function(sv, flank_del=50) {
  
  # only get BND records
  tmp=sv %>%
    filter(grepl("BND", INFO)) %>%
    mutate(chr2=paste0('chr',gsub('.*chr(.*):.*','\\1', ALT)),
           pos2=gsub('.*chr.*:(\\d+).*','\\1', ALT),
           pos2=as.integer(pos2))
  
  # get BND information
  bnd_tmp=tmp %>%
    mutate(mateid=gsub(".*MATEID=([^;]+);.*", "\\1", INFO))%>%
    rowwise()%>%
    mutate(left=which(CHROM==tmp$CHROM & abs(POS-tmp$POS)<400 & chr2==tmp$chr2)[1]*2-1, # check left breakpoint to group
           right=which(chr2==tmp$chr2 & abs(pos2-tmp$pos2)<400 & CHROM==tmp$CHROM)[1]*2-1) %>% # check right breakpoint to group
    group_by(left)%>%
    mutate(grp_right=sum(right))%>% # group based on left breakpoint
    group_by(right)%>%
    mutate(grp_left=sum(left))%>% # group based on right breakpoint
    group_by(grp_left, grp_right)%>% # with left and right group, the left and right breakpoints are matched, then find stats
    mutate(len1=max(POS)-min(POS),
           len2=max(pos2)-min(pos2),
           donor=unique(CHROM[which.max(len1)]),
           receiver=unique(chr2[which.max(len2)]),
           len_dif=abs(len1-len2)) %>%
    rowwise()%>%
    mutate(match_pos_id=which(abs(pos2-tmp$POS)<50)[1], 
           match_chrom=(chr2==tmp$CHROM[match_pos_id])) %>% # find which two records are mates: if they are, pos2 should have a match in POS and same chromosome
    ungroup()%>%
    mutate(match_pos=ifelse(match_chrom==TRUE, grp_left[match_pos_id], NA),
           mate=as.character(match_pos+grp_right))%>% #label mate recods
    group_by(len_dif)%>%
    mutate(class=ifelse(max(len_dif)>50, 'utl','rtl')) %>%
    filter(!is.na(mate))
  
  # Identify DEL overlapping BND
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
          abs(POS-bnd_tmp$pos2)<flank_del
      )
    )%>%
    filter(overlap) %>%
    # assign trans id to associated del
    mutate(
      BNDid = bnd_tmp$ID[which(
        CHROM == bnd_tmp$chr2 &
          abs(POS-bnd_tmp$pos2)<flank_del
      )]
    ) %>%
    ungroup()
  
  # Annotate BND classes
  bnd=bnd_tmp %>%
    # use del to classify cut-paste, and rest to copy paste
    mutate(
      class = case_when(
        ID %in% del$BNDid ~ 'cuttl',
        class == 'utl'      ~ 'coptl',
        TRUE                 ~ class
      ))%>%
    group_by(len_dif)%>%
    mutate(class=ifelse(any(grepl('cuttl', class)), 'cuttl',class)) %>%
    group_by(class, mate, CHROM) %>%
    mutate(nPOS=min(POS),
           mPOS=max(POS))%>%
    group_by(mate)%>%
    mutate(n=n())%>%
    filter(n>1)%>%
    mutate(type=ifelse(class=='rtl' & n<4, 'out','stay'))%>%
    filter(type=='stay')%>%
    ungroup()
  
  return(list(bnd, del))
}

# # 2b. Identify DEL overlapping BND
# process_del <- function(sv, bnd_tmp, flank_del = 50) {
# 
#   # isolate deletion to find cut-paste trans
#   # check if a deletion overlaps with a trans
#   del=sv %>%
#     filter(str_detect(INFO, 'DEL'))%>%
#     mutate(chr2=CHROM,
#            pos2=gsub('.*END=(.*);CIPOS.*','\\1',INFO),
#            pos2=as.integer(pos2))%>%
#     rowwise() %>%
#     mutate(
#       #check any del is associated with trans (cut-paste)
#       overlap = any(
#         CHROM == bnd_tmp$chr2 &
#           abs(POS-bnd_tmp$pos2)<flank_del
#       )
#     )%>%
#     filter(overlap) %>%
#     # assign trans id to associated del
#     mutate(
#       BNDid = bnd_tmp$ID[which(
#         CHROM == bnd_tmp$chr2 &
#           abs(POS-bnd_tmp$pos2)<flank_del
#       )]
#     ) %>%
#     ungroup()
#   return(del)
# }
# 
# # 2c. Annotate BND classes
# annotate_bnd <- function(bnd_tmp, del) {
# 
#   # assign cut-paste, copy-paste, and reciprocal
#   bnd=bnd_tmp %>%
#     # use del to classify cut-paste, and rest to copy paste
#     mutate(
#       class = case_when(
#         ID %in% del$BNDid ~ 'cuttl',
#         class == 'utl'      ~ 'coptl',
#         TRUE                 ~ class
#       ))%>%
#     group_by(len_dif)%>%
#     mutate(class=ifelse(any(grepl('cuttl', class)), 'cuttl',class)) %>%
#     group_by(class, mate, CHROM) %>%
#     mutate(nPOS=min(POS),
#            mPOS=max(POS))%>%
#     group_by(mate)%>%
#     mutate(n=n())%>%
#     filter(n>1)%>%
#     mutate(type=ifelse(class=='rtl' & n<4, 'out','stay'))%>%
#     filter(type=='stay')%>%
#     ungroup()
# 
#   return(bnd)
# }