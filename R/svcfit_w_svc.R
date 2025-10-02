#' Run SVCFit with reads from SVclone
#'
#' @param svc an object of class 'data frame'. This object stores the output from `analyze_svclone`.
#' @param final_truth an object of class 'data frame'. This object stores the output from `calc_svcf`.
#' @param thresh an object of class 'numeric'. This object is strictly between 0 and 1, and describes the
#' threshold for SVCFit to decide whether an SV occurs before CNV or after.
#' @param sampan object of class 'character'. This object describes the sample name.
#' @param exper an object of class 'character'. This object describes the experiment number.
#'
#' @returns sth
#' @export
#'
svcfit_w_svc <-function(svc, final_truth, thresh, samp, exper){

  w_svc_read = svc %>%
    select(ID, ref_counts1, ref_counts2, var_counts1, var_counts2)%>%
    left_join(final_truth) %>%
    mutate(sv_ref=ref_counts1+ref_counts2,
           sv_alt=var_counts1)%>%
    rowwise()%>%
    mutate(raw_svcf=ifelse(classification=="DUP", major*sv_alt/(sv_alt+sv_ref), 2*sv_alt/(sv_alt+sv_ref)),
           raw_svcf=ifelse(zygosity=='hom', raw_svcf/2, raw_svcf),
           s1=(2*sv_alt-sv_ref*(ASCN-1))/(sv_alt+sv_ref),
           s2=((sv_alt+sv_alt*ASCN)/(sv_alt+sv_ref))%%2,
           ss1=ifelse(zygosity=="hom", 0.5*s1, s1),
           ss2=ifelse(zygosity=="hom", 0.5*s2, s2),
           final_svcf=ifelse(cn_type=="DEL" | ss1<=thresh, ss2,ss1 ),
           final_svcf=ifelse(classification=="COPBND" & cn_type=='DUP', round(cncf, 2), final_svcf),
           ## if del or dup is detected in facet, directly use cncf to avoid unnecessary modification when cnv doesn't overlap with them
           final_svcf=ifelse(cn_type=='DEL' & classification=='DEL' & sv_phase==cnv_phase, round(cncf, 2), round(final_svcf,2)),
           final_svcf=ifelse(cn_type=='DUP' & classification=='DUP' & sv_phase==cnv_phase, round(cncf, 2), round(final_svcf,2)),
           expmt=exper,
           sample=samp)%>%
    ungroup()

  return(w_svc_read)
}
