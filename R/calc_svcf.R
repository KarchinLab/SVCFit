#' Calculate SVCF
#'
#' @param anno_sv_cnv an object of class 'data frame'. This object stores the output from `annotate_cnv`.
#' @param sv_info an object of class 'data frame'. This object stores the output from `parse_sv_info`.
#' @param thresh an object of class 'numeric'. This object is strictly between 0 and 1, and describes the
#' threshold for SVCFit to decide whether an SV occurs before CNV or after.
#' @param samp an object of class 'character'. This object describes the sample name.
#' @param exper an object of class 'character'. This object describes the experiment number.
#'
#' @returns sth
#' @export
#'
calc_svcf <- function(anno_sv_cnv, sv_info, thresh=0.1, samp, exper){
  final <- inner_join(anno_sv_cnv, sv_info)%>%
    filter(!is.na(ASCN))%>%
    ##ss2 is cnv and sv diff phase, ss1 is same phase
    mutate(raw_svcf=ifelse(classification=="DUP", major*sv_alt/(sv_alt+sv_ref), 2*sv_alt/(sv_alt+sv_ref)),
           raw_svcf=ifelse(zygosity=='hom', raw_svcf/2, raw_svcf),
           s1=(2*sv_alt-sv_ref*(ASCN-1))/(sv_alt+sv_ref),
           s2=((sv_alt+sv_alt*ASCN)/(sv_alt+sv_ref))%%2,
           ss1=ifelse(zygosity=="hom", 0.5*s1, s1),
           ss2=ifelse(zygosity=="hom", 0.5*s2, s2),
           ## decide order (only consider del happened before sv, so cn_type==del always ss2)
           final_svcf=ifelse(cn_type=="DEL" | ss1<=thresh, ss2,ss1 ),
           final_svcf=ifelse(classification=="COPBND" & cn_type=='DUP', round(cncf, 2), final_svcf),
           ## if del or dup is detected in facet, directly use cncf to avoid unnecessary modification when cnv doesn't overlap with them
           final_svcf=ifelse(cn_type=='DEL' & classification=='DEL' , round(cncf, 2), round(final_svcf,2)),
           final_svcf=ifelse(cn_type=='DUP' & classification=='DUP' , round(cncf, 2), round(final_svcf,2)),
           expmt=exper,
           sample=samp,
           tmp_id=ifelse(grepl('BND', ID), gsub("((:[^:\\D:]){6}$)","", ID),gsub("((:[^:\\D:]){5}$)","", ID)))%>%
    group_by(mate)%>%
    mutate(final_svcf=mean(final_svcf),
           classification=ifelse(grepl('BND', classification), 'BND', classification))%>%
    ungroup()
  
  return(final)
}
