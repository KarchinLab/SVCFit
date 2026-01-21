#' Load and process data for clustering
#'
#' Reads paired sample IDs and purity files, calculates Cancer Cell Fraction (CCF) 
#' for Structural Variants (SVs), and reshapes the data (wide format) for longitudinal clustering.
#' It handles shared and non-shared variants between pre- and on-treatment samples.
#'
#' @param pair_path Character. Path to a file containing paired sample IDs (columns: pre_BAT, on_BAT).
#' @param pur_path Character. Path to a file containing purity estimates for each sample.
#'
#' @return A data frame containing processed CCF values for pre- and on-treatment samples 
#' (`pre_BAT`, `on_BAT`), ready for input into clustering functions.
#' @export
pre_process_cluster <- function(pair_path, pur_path){
  sample_pair <- read.delim(pair_path, header=FALSE)
  colnames(sample_pair)=c("pre_BAT","on_BAT")
  samp_pur <- read.delim(pur_path)
  combine_data = lapply(1:nrow(sample_pair), function(x) read_data(sample_pair, x, samp_pur))%>%
    do.call(rbind,.)%>%
    ungroup()%>%
    filter(pair!=9)%>%
    mutate(ccf=final_svcf/purity,
           ccf=ifelse(ccf>1, 1, ccf),
           final_svcf=ifelse(final_svcf>1, 1, final_svcf))%>%
    #filter(!sample_ID %in% responder_id)%>%
    mutate(pair=as.character(pair),
           sv_len=mPOS-nPOS)%>%
    group_by(pair, sid)%>%
    mutate(POS=ifelse(!is.na(shared),mean(POS), POS),
           END=ifelse(!is.na(shared),mean(END), END),
           sv_len=mean(sv_len))%>%
    ungroup()%>%
    filter(sv_len>50)
  
  shared=combine_data %>%
    filter(!is.na(shared))%>%
    select(CHROM, POS, END, ccf, pair, classification, stage, sid, sv_len)%>%
    group_by(CHROM, POS, stage)%>%
    mutate(ccf=round(mean(ccf), digits=2))%>%
    distinct(.keep_all = T)%>%
    ungroup() %>%
    filter(!is.na(sid))%>%
    pivot_wider(names_from = stage, values_from = ccf)%>%
    rowwise()%>%
    mutate(new_pre = mean(pre_BAT[[1]]),
           new_on = mean(on_BAT[[1]]),
           pre_BAT=new_pre,
           on_BAT=new_on)%>%
    ungroup()%>%
    select(-c('new_pre','new_on'))
  
  non_shared=combine_data %>% 
    filter(is.na(shared)) %>%
    select(CHROM, POS, END, ccf, pair, classification, stage,sid, sv_len)%>%
    group_by(CHROM, POS, stage)%>%
    mutate(ccf=round(mean(ccf), digits=2))%>%
    distinct(.keep_all = T)%>%
    ungroup() %>%
    pivot_wider(names_from = stage, values_from = ccf) %>%
    rowwise()%>%
    mutate(pre_BAT = ifelse(is.na(pre_BAT), 0, pre_BAT),
           on_BAT = ifelse(is.na(on_BAT), 0, on_BAT))%>%
    ungroup()
  
  new_dat <- rbind(shared, non_shared) %>%
    filter((pre_BAT+on_BAT)!=0) %>%
    mutate(
      pre_BAT = ifelse(pre_BAT < 0.1, 0, pre_BAT),
      on_BAT = ifelse(on_BAT < 0.1, 0, on_BAT)
    )
  return(new_dat)
}

read_data <- function(pair, row, pur_file){
  pre_samp=pair$pre_BAT[row]
  on_samp=pair$on_BAT[row]
  pre_pur=pur_file$purity[pur_file$sample==pre_samp]
  on_pur=pur_file$purity[pur_file$sample==on_samp]
  pre <- read.delim(paste0(mendeley_dir, '/COMBAT/SVCFit_output/', pre_samp, '.bed')) %>%
    mutate(sample_ID=pre_samp,
           pair=row,
           stage='pre_BAT',
           purity=as.numeric(pre_pur))
  on <- read.delim(paste0(mendeley_dir, '/COMBAT/SVCFit_output/', on_samp, '.bed')) %>%
    mutate(sample_ID=on_samp,
           pair=row,
           stage='on_BAT',
           purity=as.numeric(on_pur))%>%
    rowwise()%>%
    mutate(row=ifelse(any(CHROM==pre$CHROM & abs(POS-pre$POS)<100 & abs(END-pre$END)<100 & classification==pre$classification),which(CHROM==pre$CHROM & abs(POS-pre$POS)<100 & abs(END-pre$END)<100 & classification==pre$classification), 0),
           shared=ifelse(row==0, NA, pre$ID[row]))
  shared_id=on$shared[!is.na(on$shared)]
  pre=pre %>%
    rowwise()%>%
    mutate(row=ifelse(any(CHROM==on$CHROM & abs(POS-on$POS)<100 & abs(END-on$END)<100 & classification==on$classification),which(CHROM==on$CHROM & abs(POS-on$POS)<100 & abs(END-on$END)<100 & classification==on$classification), 0),
           shared=ifelse(row==0, NA, on$ID[row]),
           sid=ifelse(ID %in% shared_id, which(ID==shared_id), NA))%>%
    ungroup()
  on <- on %>%
    mutate(sid=ifelse(shared %in% shared_id, which(shared==shared_id), NA))%>%
    ungroup()
  
  out = rbind(pre, on) %>% ungroup()
  
  return(out)
}


##clustering util

dp_gmm_convergence <- function(Z, Kmax = 10, n_steps = 50, random_state = 0L, concentration) {
  
  bgmm <- sk$mixture$BayesianGaussianMixture(
    n_components = as.integer(Kmax),
    covariance_type = "full",
    weight_concentration_prior_type = "dirichlet_process",
    weight_concentration_prior = concentration,
    init_params = "random",
    n_init = as.integer(1),      # important for warm_start
    max_iter = as.integer(1),    # one EM step per fit()
    reg_covar = 1e-6,
    random_state = as.integer(random_state),
    warm_start = TRUE
  )
  
  lower_bounds <- numeric(n_steps)
  
  for (i in seq_len(n_steps)) {
    bgmm$fit(Z)
    lower_bounds[i] <- bgmm$lower_bound_
  }
  
  tibble(
    iter        = seq_len(n_steps),
    lower_bound = lower_bounds
  )
}

logit_eps <- function(p, n = 200) {
  q <- (p * (n - 1) + 0.5) / n   # principled epsilon for 0/1
  log(q / (1 - q))
}

#' Run Dirichlet Process Gaussian Mixture Model (DP-GMM) clustering
#'
#' Fits a Bayesian Gaussian Mixture Model (via reticulate and sklearn) to a specific pair of samples.
#' It models the distribution of CCF changes (or absolute values) to identify mutation clusters.
#' Includes convergence checks and plotting.
#'
#' @param input Data frame. The processed output from `pre_process_cluster`.
#' @param pair_num Character or Numeric. The identifier/row number for the specific sample pair to cluster.
#' @param Kmax Numeric. The maximum number of clusters (components) to allow in the Dirichlet Process. Default is 10.
#' @param n_steps Numeric. The number of steps for the convergence check run. Default is 100.
#' @param thr_min_w Numeric. The minimum weight threshold for a component to be considered active. Default is 0.01.
#' @param random_state Integer. Seed for the random number generator in Python/sklearn for reproducibility.
#' @param concentration Numeric. The Dirichlet concentration prior (gamma). Higher values encourage more active components; lower values encourage sparsity. Default is 1.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{clusters}: Data frame with cluster assignments and probabilities.
#'   \item \code{model}: The fitted sklearn BayesianGaussianMixture object.
#'   \item \code{plot_conv}: A ggplot object showing the convergence (ELBO) trace.
#'   \item \code{plot_clust}: A ggplot object showing the scatter plot of the clusters.
#' }
#' @export
run_dp_gmm_pair <- function(input, pair_num,
                            Kmax      = 10,
                            n_steps   = 100,
                            thr_min_w = 0.01,
                            random_state = 0L,
                            concentration=1) {
  # 1) Subset and prepare data
  X <- input %>%
    filter(pair == pair_num) %>%
    distinct(pre_BAT, on_BAT, .keep_all = TRUE)
  
  if (nrow(X) == 0) {
    stop("No rows found for this pair_num.")
  }
  
  Z <- cbind(
    (X$on_BAT-X$pre_BAT),
    (X$pre_BAT)
  )
  
  # Z <- cbind(logit_eps(X$pre_BAT)/logit_eps(X$on_BAT),
  #            X$pre_BAT)
  
  ## 2) Convergence run (short, for plotting)
  conv_df <- dp_gmm_convergence(
    Z,
    Kmax = Kmax,
    n_steps = n_steps,
    random_state = random_state,
    concentration
  )
  
  p_conv <- ggplot(conv_df, aes(x = iter, y = lower_bound)) +
    geom_line() +
    geom_point() +
    labs(
      x = "EM iteration",
      y = "Variational lower bound",
      title = paste("DP-GMM convergence (pair", pair_num, ")")
    ) +
    theme_classic()
  
  # Optional: change in lower bound
  conv_df_delta <- conv_df %>%
    arrange(iter) %>%
    mutate(delta = c(NA, diff(lower_bound)))
  
  p_delta <- ggplot(conv_df_delta, aes(x = iter, y = delta)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line() +
    geom_point() +
    labs(
      x = "EM iteration",
      y = "Change in lower bound",
      title = paste("Change in lower bound per iteration (pair", pair_num, ")")
    ) +
    theme_classic()
  
  ## 3) Final model fit (more iterations, multiple inits) for clustering
  bgmm_final <- sk$mixture$BayesianGaussianMixture(
    n_components = as.integer(Kmax),
    covariance_type = "full",
    weight_concentration_prior_type = "dirichlet_process",
    weight_concentration_prior = concentration,
    init_params = "random",
    n_init = as.integer(5),
    max_iter = as.integer(2000),
    reg_covar = 1e-6,
    random_state = as.integer(random_state)
  )$fit(Z)
  
  cat("converged:", bgmm_final$converged_, "\n")
  cat("n_iter:", bgmm_final$n_iter_, "\n")
  cat("final lower_bound:", bgmm_final$lower_bound_, "\n")
  
  ## 4) Responsibilities and active components
  resp    <- bgmm_final$predict_proba(Z)      # N x Kmax
  weights <- as.numeric(bgmm_final$weights_)  # mixture weights
  
  thr_w <- max(1 / nrow(Z), thr_min_w)        # threshold for "active" comps
  active <- which(weights > thr_w)
  
  if (length(active) == 0) {
    warning("No active components above weight threshold; using all components.")
    active <- seq_along(weights)
  }
  
  resp_active <- resp[, active, drop = FALSE]
  hard_idx_local <- max.col(resp_active, ties.method = "first")
  hard_comp     <- active[hard_idx_local]
  hard_prob     <- apply(resp_active, 1, max)
  
  ## 5) Build clustering data frame
  df_clusters <- tibble(
    pre_ccf  = X$pre_BAT,
    post_ccf = X$on_BAT,
    cluster  = factor(match(hard_comp, active)),  # 1..#active
    prob     = pmin(pmax(hard_prob, 0), 1)        # clamp to [0,1]
  ) %>%
    group_by(cluster) %>%
    mutate(
      pre_center  = mean(pre_ccf),
      post_center = mean(post_ccf)
    ) %>%
    ungroup() %>%
    mutate(pair = pair_num)
  
  ## 6) Simple cluster scatter plot
  p_clusters <- ggplot(df_clusters,
                       aes(x = pre_ccf, y = post_ccf,
                           color = cluster, alpha = prob)) +
    geom_point() +
    labs(
      x = "pre_BAT",
      y = "on_BAT",
      title = paste("DP-GMM clusters (pair", pair_num, ")"),
      alpha = "cluster\nprob."
    ) +
    theme_classic()
  
  ## 8) Return objects for further use
  list(
    clusters   = df_clusters,
    conv_trace = conv_df,
    conv_delta = conv_df_delta,
    model      = bgmm_final,
    plot_conv  = p_conv,
    plot_delta = p_delta,
    plot_clust = p_clusters
  )
}

## post process util
#' Enforce minimum size for clusters
#'
#' Checks the size of generated clusters. If a cluster has fewer than \code{min_n} observations,
#' it is merged into the nearest valid cluster based on the Euclidean distance of their centers.
#'
#' @param df Data frame. The output containing initial cluster assignments (usually from `run_dp_gmm_pair` or merged results).
#' @param min_n Integer. The minimum number of mutations required to keep a cluster distinct. Default is 5.
#'
#' @return A data frame with updated cluster assignments in a column named \code{ncluster}.
#' @export
enforce_min_cluster_size <- function(df, min_n = 5) {
  df %>%
    mutate(cluster = as.integer(cluster)) %>%
    group_by(pair, cluster) %>%
    mutate(n_in_clust = dplyr::n()) %>%
    ungroup() %>%
    group_by(pair) %>%
    # compute one row per (pair, cluster) with centers and sizes
    do({
      cl <- .
      centers <- cl %>%
        group_by(cluster) %>%
        summarise(n = dplyr::n(),
                  pre_c = mean(pre_ccf),
                  post_c = mean(post_ccf),
                  .groups = "drop")
      
      # clusters to keep (size >= min_n); if none qualify, keep the largest as seed
      keepers <- centers %>% filter(n >= min_n)
      if (nrow(keepers) == 0) {
        keepers <- centers %>% slice_max(n, n = 1)
      }
      smalls <- centers %>%
        filter(!(cluster %in% keepers$cluster))
      
      # map small clusters → nearest keeper by center distance
      if (nrow(smalls) > 0) {
        assign_map <- smalls %>%
          rowwise() %>%
          mutate(closest_keeper = keepers$cluster[which.min(
            (pre_c - keepers$pre_c)^2 + (post_c - keepers$post_c)^2
          )]) %>%
          ungroup() %>%
          select(cluster, closest_keeper)
        
        cl <- cl %>%
          left_join(assign_map, by = "cluster") %>%
          mutate(cluster_merged = ifelse(is.na(closest_keeper), cluster, closest_keeper)) %>%
          select(-closest_keeper)
      } else {
        cl <- cl %>% mutate(cluster_merged = cluster)
      }
      cl
    }) %>%
    ungroup() %>%
    # re-compute centers after merging
    group_by(pair, cluster_merged) %>%
    mutate(
      n_after     = dplyr::n(),
      pre_center  = mean(pre_ccf),
      post_center = mean(post_ccf)
    ) %>%
    ungroup() %>%
    mutate(ncluster = factor(cluster_merged)) %>%
    select(-n_in_clust, -cluster_merged)
}

#' Identify geometrically close clusters for merging
#'
#' Calculates the Euclidean distance between cluster centers. If the euclidean distance between the centers of two clusters
#' is less than 0.2, they are flagged to be merged.
#'
#' @param input Data frame. A data frame containing `pre_center`, `post_center`, and `cluster` columns.
#' @param pair_num Character or Numeric. The identifier for the specific sample pair being analyzed.
#'
#' @return A data frame containing pairs of cluster IDs (`row`, `col`) that are close enough to be merged.
#' @export
merge_cluster <- function(input, pair_num){
  tmp = input %>% 
    filter(pair==pair_num)%>%
    select(pre_center, post_center, cluster)%>%
    arrange(cluster)%>%
    distinct(cluster, .keep_all = T)%>%
    ungroup()%>%
    select(-cluster)%>%
    as.matrix()
  d=as.matrix(dist(tmp))
  p=which(d<0.2 & d>0, arr.ind = T)%>%
    as.data.frame() %>%
    rowwise() %>%
    mutate(key = paste(sort(c_across(everything())), collapse = "_")) %>%
    ungroup() %>%
    distinct(key, .keep_all = TRUE) %>%
    select(-key)%>%
    mutate(pair=pair_num)
  return(p)
}

## cluster function
#' Perform full clustering pipeline on paired data
#'
#' This wrapper function executes the entire clustering workflow:
#' 1. Loads and processes data (`pre_process_cluster`).
#' 2. Initializes the Python environment (`reticulate`).
#' 3. Runs DP-GMM clustering on specified pairs.
#' 4. Merges geometrically close clusters.
#' 5. Enforces minimum cluster sizes.
#' 6. Maps cluster assignments back to original Structural Variant (SV) data.
#'
#' @param pair_path Character. Path to file containing paired sample IDs.
#' @param pur_path Character. Path to file containing purity for each sample.
#' @param Kmax Numeric. Maximum k for DP-GMM clustering. Default is 10.
#' @param n_steps Numeric. Maximum steps for iteration. Default is 100.
#' @param thr_min_w Numeric. Minimum weight for active components. Default is 0.01.
#' @param random_state Integer. Random seed for reproducibility. Default is 0.
#' @param concentration Numeric. Dirichlet concentration prior. Default is 1.
#' @param min_n Integer. Minimum number of mutations per cluster. Default is 5.
#' @param pair_num Character or Numeric. The identifier for the specific sample pair being analyzed.
#'
#' @return A list containing two elements:
#' \enumerate{
#'   \item \code{cluster_result}: Data frame with the final clustering centers and assignments.
#'   \item \code{sv2cluster}: Data frame mapping specific SVs (CHROM, POS, END) to their assigned clusters.
#' }
#' @export
cluster_data <- function(pair_path, 
                         pur_path, 
                         Kmax      = 10,
                         n_steps   = 100,
                         thr_min_w = 0.01,
                         random_state = 0L,
                         concentration=1,
                         min_n = 5,
                         pair_num=1){
  new_dat <- pre_process_cluster(pair_path, pur_path)
  use_condaenv("py3", required = TRUE)
  py_config()
  sk <- import("sklearn", delay_load = TRUE)
  
  df_plot1 <- lapply(c(1:8,10:12), function(x) run_dp_gmm_pair(new_dat, x)[[1]]) %>%
    do.call(rbind,.)
  
  df_plot=df_plot1 %>% mutate(cluster=as.integer(cluster))
  
  cc=lapply(c(1:8,10:12), function(x) merge_cluster(df_plot, x))%>%
    do.call(rbind,.)%>%
    dplyr::rename('cluster'='row')
  
  inter_df_plot = left_join(df_plot,cc)%>%
    mutate(ncluster=ifelse(is.na(col), cluster, col),
           ncluster=as.factor(ncluster))%>%
    group_by(pair, ncluster)%>%
    mutate(pre_center=mean(pre_ccf),
           post_center=mean(post_ccf))%>%
    ungroup()
  
  new_plot <- enforce_min_cluster_size(inter_df_plot, min_n = 5) 
  
  cc=new_plot %>%
    distinct(pair, ncluster)%>%
    arrange(pair, ncluster)%>%
    group_by(pair)%>%
    mutate(cluster_num=paste0('cluster', row_number()))
  
  cluster_result <- new_plot %>% left_join(cc)
  
  sv2cluster = cluster_result%>%
    mutate(pre_BAT=pre_ccf,
           on_BAT=post_ccf,
           pair=as.character(pair))%>%
    right_join(new_dat, by=c('pair','pre_BAT','on_BAT'))%>%
    mutate(sv_len=END-POS)%>%
    select(CHROM, POS, END, classification, pair, cluster_num, pre_center, post_center, pre_BAT, on_BAT, sv_len)
  
  clones = cluster_result %>%
    select(pair, cluster_num, pre_center, post_center)%>%
    group_by(pair, cluster_num)%>%
    distinct(.keep_all = T)%>%
    ungroup()%>%
    mutate(f_pre=pre_center, f_day85=post_center,
           se_pre=0.1, se_day85=0.1,
           var_pre = se_pre^2,
           var_day85 = se_day85^2,
           w_pre = 1/var_pre,
           w_day85 = 1/var_day85)%>%
    select(pair, cluster_num, f_pre, f_day85, se_pre, se_day85, var_pre, var_day85, w_pre, w_day85)%>%
    filter(pair==pair_num)%>%
    arrange(cluster_num)
  
  return(list(cluster_result, sv2cluster, clones))
}


