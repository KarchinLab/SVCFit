#' Load and process data for clustering
#'
#' Reads paired sample IDs and purity files, calculates Cancer Cell Fraction (CCF) 
#' for Structural Variants (SVs), and reshapes the data (wide format) for longitudinal clustering.
#' It handles shared and non-shared variants between pre- and on-treatment samples.
#'
#' Key fix with zero effect on clustering:
#' - Keep the original clustering workflow unchanged
#' - Still overwrite POS/END with mean breakpoint for shared events in the clustering table
#' - Preserve original raw coordinates separately for mapping/AnnotSV
#' - Add event_id so cluster labels can be merged back to original SV rows
#'
#' @param pair_path Character. Path to a file containing paired sample IDs (columns: pre_BAT, on_BAT).
#' @param pur_path Character. Path to a file containing purity estimates for each sample.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{cluster_input}: event-level table for clustering (same logic as original)
#'   \item \code{sv_raw}: original row-level SV table with raw coordinates
#' }
#' @export
pre_process_cluster <- function(pair_path, pur_path){
  sample_pair <- read.delim(pair_path, header = FALSE)
  colnames(sample_pair) <- c("pre_BAT", "on_BAT")
  samp_pur <- read.delim(pur_path)
  
  combine_data <- lapply(1:nrow(sample_pair), function(x) read_data(sample_pair, x, samp_pur)) %>%
    do.call(rbind, .) %>%
    ungroup() %>%
    filter(pair != 9) %>%
    mutate(
      ccf = final_svcf / purity,
      ccf = ifelse(ccf > 1, 1, ccf),
      final_svcf = ifelse(final_svcf > 1, 1, final_svcf),
      pair = as.character(pair),
      sv_len = mPOS - nPOS,
      POS_raw = POS,
      END_raw = END
    ) %>%
    mutate(
      event_id = ifelse(
        !is.na(sid),
        paste0("pair", pair, "_shared_", sid),
        paste0("pair", pair, "_", sample_ID, "_", ID)
      )
    ) %>%
    group_by(pair, sid) %>%
    mutate(
      POS = ifelse(!is.na(shared), mean(POS), POS),
      END = ifelse(!is.na(shared), mean(END), END),
      sv_len = mean(sv_len)
    ) %>%
    ungroup() %>%
    filter(sv_len > 50)
  
  # keep an original-row version for later mapping / AnnotSV joins
  sv_raw <- combine_data
  
  # build the clustering input exactly the same way as the original script
  shared <- combine_data %>%
    filter(!is.na(shared)) %>%
    select(event_id, CHROM, POS, END, ccf, pair, classification, stage, sid, sv_len) %>%
    group_by(CHROM, POS, stage) %>%
    mutate(ccf = round(mean(ccf), digits = 2)) %>%
    distinct(.keep_all = TRUE) %>%
    ungroup() %>%
    filter(!is.na(sid)) %>%
    pivot_wider(names_from = stage, values_from = ccf) %>%
    rowwise() %>%
    mutate(
      new_pre = mean(pre_BAT[[1]]),
      new_on = mean(on_BAT[[1]]),
      pre_BAT = new_pre,
      on_BAT = new_on
    ) %>%
    ungroup() %>%
    select(-c("new_pre", "new_on"))
  
  non_shared <- combine_data %>% 
    filter(is.na(shared)) %>%
    select(event_id, CHROM, POS, END, ccf, pair, classification, stage, sid, sv_len) %>%
    group_by(CHROM, POS, stage) %>%
    mutate(ccf = round(mean(ccf), digits = 2)) %>%
    distinct(.keep_all = TRUE) %>%
    ungroup() %>%
    pivot_wider(names_from = stage, values_from = ccf) %>%
    rowwise() %>%
    mutate(
      pre_BAT = ifelse(is.na(pre_BAT), 0, pre_BAT),
      on_BAT = ifelse(is.na(on_BAT), 0, on_BAT)
    ) %>%
    ungroup()
  
  new_dat <- rbind(shared, non_shared) %>%
    filter((pre_BAT + on_BAT) != 0) %>%
    mutate(
      pre_BAT = ifelse(pre_BAT < 0.1, 0, pre_BAT),
      on_BAT = ifelse(on_BAT < 0.1, 0, on_BAT),
      pair = as.character(pair)
    )
  
  return(list(
    cluster_input = new_dat,
    sv_raw = sv_raw
  ))
}

read_data <- function(pair, row, pur_file){
  pre_samp <- pair$pre_BAT[row]
  on_samp <- pair$on_BAT[row]
  pre_pur <- pur_file$purity[pur_file$sample == pre_samp]
  on_pur <- pur_file$purity[pur_file$sample == on_samp]
  
  pre <- read.delim(paste0(mendeley_dir, '/COMBAT/SVCFit_output/', pre_samp, '.bed')) %>%
    mutate(
      sample_ID = pre_samp,
      pair = row,
      stage = 'pre_BAT',
      purity = as.numeric(pre_pur)
    )
  
  on <- read.delim(paste0(mendeley_dir, '/COMBAT/SVCFit_output/', on_samp, '.bed')) %>%
    mutate(
      sample_ID = on_samp,
      pair = row,
      stage = 'on_BAT',
      purity = as.numeric(on_pur)
    ) %>%
    rowwise() %>%
    mutate(
      row = ifelse(
        any(CHROM == pre$CHROM & abs(POS - pre$POS) < 100 & abs(END - pre$END) < 100 & classification == pre$classification),
        which(CHROM == pre$CHROM & abs(POS - pre$POS) < 100 & abs(END - pre$END) < 100 & classification == pre$classification)[1],
        0
      ),
      shared = ifelse(row == 0, NA, pre$ID[row])
    )
  
  shared_id <- on$shared[!is.na(on$shared)]
  
  pre <- pre %>%
    rowwise() %>%
    mutate(
      row = ifelse(
        any(CHROM == on$CHROM & abs(POS - on$POS) < 100 & abs(END - on$END) < 100 & classification == on$classification),
        which(CHROM == on$CHROM & abs(POS - on$POS) < 100 & abs(END - on$END) < 100 & classification == on$classification)[1],
        0
      ),
      shared = ifelse(row == 0, NA, on$ID[row]),
      sid = ifelse(ID %in% shared_id, which(ID == shared_id), NA)
    ) %>%
    ungroup()
  
  on <- on %>%
    mutate(
      sid = ifelse(shared %in% shared_id, which(shared == shared_id), NA)
    ) %>%
    ungroup()
  
  out <- rbind(pre, on) %>% ungroup()
  return(out)
}


## clustering util

dp_gmm_convergence <- function(Z, Kmax = 10, n_steps = 50, random_state = 0L, concentration) {
  
  bgmm <- sk$mixture$BayesianGaussianMixture(
    n_components = as.integer(Kmax),
    covariance_type = "full",
    weight_concentration_prior_type = "dirichlet_process",
    weight_concentration_prior = concentration,
    init_params = "random",
    n_init = as.integer(1),
    max_iter = as.integer(1),
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
    iter = seq_len(n_steps),
    lower_bound = lower_bounds
  )
}

logit_eps <- function(p, n = 200) {
  q <- (p * (n - 1) + 0.5) / n
  log(q / (1 - q))
}

#' Run Dirichlet Process Gaussian Mixture Model (DP-GMM) clustering
#'
#' This function is intentionally kept as close as possible to the original,
#' to preserve clustering behavior.
#'
#' @export
run_dp_gmm_pair <- function(input, pair_num,
                            Kmax = 10,
                            n_steps = 100,
                            thr_min_w = 0.01,
                            random_state = 0L,
                            concentration = 1) {
  X <- input %>%
    filter(pair == pair_num) %>%
    distinct(pre_BAT, on_BAT, .keep_all = TRUE)
  
  if (nrow(X) == 0) {
    stop("No rows found for this pair_num.")
  }
  
  Z <- cbind(
    (X$on_BAT - X$pre_BAT),
    (X$pre_BAT)
  )
  
  conv_df <- dp_gmm_convergence(
    Z,
    Kmax = Kmax,
    n_steps = n_steps,
    random_state = random_state,
    concentration = concentration
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
  
  resp <- bgmm_final$predict_proba(Z)
  weights <- as.numeric(bgmm_final$weights_)
  
  thr_w <- max(1 / nrow(Z), thr_min_w)
  active <- which(weights > thr_w)
  
  if (length(active) == 0) {
    warning("No active components above weight threshold; using all components.")
    active <- seq_along(weights)
  }
  
  resp_active <- resp[, active, drop = FALSE]
  hard_idx_local <- max.col(resp_active, ties.method = "first")
  hard_comp <- active[hard_idx_local]
  hard_prob <- apply(resp_active, 1, max)
  
  df_clusters <- tibble(
    pre_ccf = X$pre_BAT,
    post_ccf = X$on_BAT,
    cluster = factor(match(hard_comp, active)),
    prob = pmin(pmax(hard_prob, 0), 1)
  ) %>%
    group_by(cluster) %>%
    mutate(
      pre_center = mean(pre_ccf),
      post_center = mean(post_ccf)
    ) %>%
    ungroup() %>%
    mutate(pair = as.character(pair_num))
  
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
  
  list(
    clusters = df_clusters,
    conv_trace = conv_df,
    conv_delta = conv_df_delta,
    model = bgmm_final,
    plot_conv = p_conv,
    plot_delta = p_delta,
    plot_clust = p_clusters
  )
}


## post process util

#' Enforce minimum size for clusters
#' @export
enforce_min_cluster_size <- function(df, min_n = 5) {
  df %>%
    mutate(cluster = as.integer(cluster)) %>%
    group_by(pair, cluster) %>%
    mutate(n_in_clust = dplyr::n()) %>%
    ungroup() %>%
    group_by(pair) %>%
    do({
      cl <- .
      centers <- cl %>%
        group_by(cluster) %>%
        summarise(
          n = dplyr::n(),
          pre_c = mean(pre_ccf),
          post_c = mean(post_ccf),
          .groups = "drop"
        )
      
      keepers <- centers %>% filter(n >= min_n)
      if (nrow(keepers) == 0) {
        keepers <- centers %>% slice_max(n, n = 1)
      }
      
      smalls <- centers %>%
        filter(!(cluster %in% keepers$cluster))
      
      if (nrow(smalls) > 0) {
        assign_map <- smalls %>%
          rowwise() %>%
          mutate(
            closest_keeper = keepers$cluster[which.min(
              (pre_c - keepers$pre_c)^2 + (post_c - keepers$post_c)^2
            )]
          ) %>%
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
    group_by(pair, cluster_merged) %>%
    mutate(
      n_after = dplyr::n(),
      pre_center = mean(pre_ccf),
      post_center = mean(post_ccf)
    ) %>%
    ungroup() %>%
    mutate(ncluster = factor(cluster_merged)) %>%
    select(-n_in_clust, -cluster_merged)
}

#' Identify geometrically close clusters for merging
#' @export
merge_cluster <- function(input, pair_num){
  tmp <- input %>%
    filter(pair == pair_num) %>%
    select(pre_center, post_center, cluster) %>%
    arrange(cluster) %>%
    distinct(cluster, .keep_all = TRUE) %>%
    ungroup() %>%
    select(-cluster) %>%
    as.matrix()
  
  if (nrow(tmp) < 2) {
    return(data.frame(row = integer(0), col = integer(0), pair = character(0)))
  }
  
  d <- as.matrix(dist(tmp))
  
  p <- which(d < 0.2 & d > 0, arr.ind = TRUE) %>%
    as.data.frame() %>%
    rowwise() %>%
    mutate(key = paste(sort(c_across(everything())), collapse = "_")) %>%
    ungroup() %>%
    distinct(key, .keep_all = TRUE) %>%
    select(-key) %>%
    mutate(pair = as.character(pair_num))
  
  return(p)
}


## cluster function

#' Perform full clustering pipeline on paired data
#'
#' This preserves the original clustering path and only changes the
#' final mapping back to original SV rows.
#'
#' @export
cluster_data <- function(pair_path, 
                         pur_path, 
                         Kmax = 10,
                         n_steps = 100,
                         thr_min_w = 0.01,
                         random_state = 0L,
                         concentration = 1,
                         min_n = 5,
                         pair_num = 1){
  
  prep <- pre_process_cluster(pair_path, pur_path)
  new_dat <- prep$cluster_input
  sv_raw <- prep$sv_raw
  
  #use_condaenv("py3", required = TRUE)
  py_config()
  sk <<- import("sklearn", delay_load = TRUE)
  
  # keep original pair handling to preserve behavior
  df_plot1 <- lapply(c(1:8, 10:12), function(x) run_dp_gmm_pair(new_dat, x)[[1]]) %>%
    do.call(rbind, .)
  
  df_plot <- df_plot1 %>%
    mutate(
      cluster = as.integer(cluster),
      pair = as.character(pair)
    )
  
  cc_list <- lapply(c(1:8, 10:12), function(x) merge_cluster(df_plot, x))
  cc <- do.call(rbind, cc_list)
  
  if (nrow(cc) > 0) {
    cc <- cc %>% dplyr::rename(cluster = row)
    
    inter_df_plot <- left_join(df_plot, cc, by = c("cluster", "pair")) %>%
      mutate(
        ncluster = ifelse(is.na(col), cluster, col),
        ncluster = as.factor(ncluster)
      ) %>%
      group_by(pair, ncluster) %>%
      mutate(
        pre_center = mean(pre_ccf),
        post_center = mean(post_ccf)
      ) %>%
      ungroup()
  } else {
    inter_df_plot <- df_plot %>%
      mutate(ncluster = as.factor(cluster)) %>%
      group_by(pair, ncluster) %>%
      mutate(
        pre_center = mean(pre_ccf),
        post_center = mean(post_ccf)
      ) %>%
      ungroup()
  }
  
  new_plot <- enforce_min_cluster_size(inter_df_plot, min_n = min_n)
  
  cc <- new_plot %>%
    distinct(pair, ncluster) %>%
    arrange(pair, ncluster) %>%
    group_by(pair) %>%
    mutate(cluster_num = paste0("cluster", row_number())) %>%
    ungroup()
  
  cluster_result <- new_plot %>%
    left_join(cc, by = c("pair", "ncluster"))
  
  # -------------------------------------------------------------------------
  # Map the cluster assignment from each UNIQUE (pair, pre_BAT, on_BAT) point
  # back to ALL original events sharing that point.
  #
  # This preserves the clustering result because clustering is still done on
  # the deduplicated coordinates only. We are only propagating the resulting
  # label back to all original rows with the same coordinates.
  # -------------------------------------------------------------------------
  cluster_point <- cluster_result %>%
    mutate(pair = as.character(pair)) %>%
    select(pair, pre_ccf, post_ccf, cluster_num, pre_center, post_center) %>%
    distinct()
  
  cluster_event <- new_dat %>%
    mutate(pair = as.character(pair)) %>%
    left_join(
      cluster_point,
      by = c("pair", "pre_BAT" = "pre_ccf", "on_BAT" = "post_ccf")
    ) %>%
    select(
      event_id,
      pair,
      cluster_num,
      pre_center,
      post_center,
      pre_ccf = pre_BAT,
      post_ccf = on_BAT
    ) %>%
    distinct()
  
  sv2cluster <- sv_raw %>%
    mutate(pair = as.character(pair)) %>%
    left_join(
      cluster_event,
      by = c("event_id", "pair")
    ) %>%
    mutate(
      sv_len = END_raw - POS_raw,
      pre_BAT = pre_ccf,
      on_BAT = post_ccf,
      POS_cluster = POS,
      END_cluster = END
    ) %>%
    select(
      event_id,
      sample_ID,
      stage,
      CHROM,
      POS = POS_raw,
      END = END_raw,
      POS_cluster,
      END_cluster,
      classification,
      pair,
      cluster_num,
      pre_center,
      post_center,
      pre_BAT,
      on_BAT,
      sv_len,
      shared,
      sid
    )
  
  clones <- cluster_result %>%
    select(pair, cluster_num, pre_center, post_center) %>%
    group_by(pair, cluster_num) %>%
    distinct(.keep_all = TRUE) %>%
    ungroup() %>%
    mutate(
      f_pre = pre_center,
      f_day85 = post_center,
      se_pre = 0.1,
      se_day85 = 0.1,
      var_pre = se_pre^2,
      var_day85 = se_day85^2,
      w_pre = 1 / var_pre,
      w_day85 = 1 / var_day85
    ) %>%
    select(pair, cluster_num, f_pre, f_day85, se_pre, se_day85,
           var_pre, var_day85, w_pre, w_day85) %>%
    filter(pair == as.character(pair_num)) %>%
    arrange(cluster_num)
  
  return(list(cluster_result, sv2cluster, clones))
}
