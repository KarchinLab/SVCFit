#' Process output from SVclone
#'
#' @param samp an object of class 'character'. This object describes sample name.
#' @param exper an object of class 'character'. This object describes experiment number.
#' @param final_truth an object of class 'data frame'. This object stores the output from SVCFit.
#' @param sv_info an object of class 'data frame'. This object stores output from `parse_sv_info`
#' @param tolerance an object of class 'integer'. This object describes the limit of
#' base pair difference between svclone and svcfit output to be considered the same event
#'
#' @returns sth
#' @export
#'
analyze_svclone <- function(samp, exper, final_truth, sv_info, tolerance=100){
  load(paste0("~/Documents/Karchin_lab/visor/m_benchmark/svclone/",samp,"/",exper,"/ccube_out/",exper,"_ccube_sv_results.RData"))
  pur=mean(doubleBreakPtsRes$ssm$purity)
  svc=doubleBreakPtsRes$ssm%>%
    mutate(first=gsub("^(.*)_.*","\\1",mutation_id),
           second=gsub("^.*_(.*)","\\1",mutation_id),
           CHROM=gsub("(.*):\\d+.*","\\1",first),
           POS=as.numeric(gsub(".*:(\\d+).*","\\1",first)),
           chr2=gsub("(.*):\\d+.*","\\1",second),
           pos2=as.numeric(gsub(".*:(\\d+).*","\\1",second)),
           ## here, try to select the same mutation present in SVCFit
           #type = any(which(abs(POS-in_dat$POS)<tolerance & abs(pos2-in_dat$END)<tolerance)),
           type = any(which(abs(POS-sv_info$POS)<tolerance)),
           row=ifelse(type == TRUE,which(abs(POS-sv_info$POS)<tolerance), 0),
           ID=ifelse(row==0,"None",sv_info$ID[row]),
           tmpid=gsub("((:[^:\\D:]){4}$)","", ID),
           svcf=(ccube_ccf1+ccube_ccf2)/2*pur,
           sample=samp,
           expmt=exper)%>%
    filter(any(grepl(tmpid, final_truth$ID)))%>% ## here, force the svclone SV to match SVCFit
    mutate(ID=final_truth$ID[which(grepl(tmpid, final_truth$ID))[1]],
           zygosity=final_truth$zygosity[which(ID==final_truth$ID)],
           classification=final_truth$classification[which(ID==final_truth$ID)],
           clone=final_truth$clone[which(ID==final_truth$ID)])%>%
    distinct(ID, .keep_all=T)
  return(svc)
}
