# Script to understand behavior of online algs as prop of non-nulls changes
# Lathan Liou

task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))

library(onlineFDR)
library(tidyverse)
library(edgeR)
library(compcodeR)
library(qvalue)
source("simulation_fct.R")

# Set Params ---- 
# rounding is CRUCIAL to avoid bizarre floating point errors
# super batches will have a higher pi followed by regular batches with a lower pi
pi_list <- list("0.1" = c(rep(0.3, 10), rep(0.05, 40)),
                "0.2" = c(rep(0.3, 10), rep(0.175, 40)),
                "0.3" = c(rep(0.5, 10), rep(0.25, 40)),
                "0.4" = c(rep(0.5, 10), rep(0.375, 40)),
                "0.5" = c(rep(0.5, 10), rep(0.5, 40)))

my_fc <- 1.5
my_alpha <- 0.05
n_samples <- 5
n_batches <- 50
n_in_batch <- 10000
output_dir <- "change_pi_sim/"

# Generate DE pvalues from sim rnaseq data ----

all_rnaseq_de_res <- list()

for(i in seq_along(1:length(pi_list))) {
  # generate batches of rnaseq data (each batch = 1 experimental run)
  # run DE analysis on each batch
  
  #check that batch layout equals n_batches
  
  if(!near(n_batches, length(pi_list[[i]]))) stop("Check n_batches and pi_list")
  
  rnaseq_de_res <- list()
  
  for(j in seq_along(1:length(pi_list[[i]]))) {

    rnaseq_data <- generateSyntheticData(dataset = paste0("sim", j, "-", names(pi_list)[i]),
                                         n.vars = n_in_batch,
                                         n.diffexp = pi_list[[i]][j]*n_in_batch,
                                         samples.per.cond = n_samples,
                                         repl.id = (j*15+i),
                                         minfact = 0.9,
                                         maxfact = 1.1,
                                         fraction.upregulated = 0.5,
                                         effect.size = 2^my_fc)
    rnaseq_de_res[[j]] <- run_limma_de(n_samples,
                                       rnaseq_data@count.matrix) %>%
      mutate(batch = j) %>%
      rename(id = gene,
             pval = p.value) %>%
      left_join(rnaseq_data@variable.annotations %>% rownames_to_column(var = "id"),
                by = "id") %>%
      rename(truth = differential.expression)
  }
  
  # bind all lists into a giant dataframe for further FDR control
  all_rnaseq_de_res[[i]] <- bind_rows(rnaseq_de_res)
}

names(all_rnaseq_de_res) <- names(pi_list)

# Apply FDR control ----
# Compare uncorrected (0.025, 0.05, 0.1), naive offline (at time of most recent batch),
# pooled offline (taking in all datasets), batch algs

apply_all_fdr_procedures <- function(data) {
  #fdr procedures on each data element (defined by pi) across family of batches
  n_batches <- length(unique(data$batch))
  
  #perform BH on each batch
  offline_bh_res <- data %>%
    group_by(batch) %>%
    mutate(adj_pval = p.adjust(pval, method = "BH"),
           R_off_bh = as.numeric(adj_pval < my_alpha)) %>%
    ungroup()
  
  #perform Storey BH on each batch
  offline_stbh_res <- data %>%
    group_by(batch) %>%
    mutate(adj_pval = qvalue(pval)$qvalues,
           R_off_stbh = as.numeric(adj_pval < my_alpha)) %>%
    ungroup()
  
  #all batches considered at once
  pooled_offline_bh_res <- data %>%
    mutate(adj_pval = p.adjust(pval, method = "BH"),
           R_pool_bh = as.numeric(adj_pval < my_alpha))
  
  pooled_offline_stbh_res <- data %>%
    mutate(adj_pval = qvalue(pval)$qvalues,
           R_pool_stbh = as.numeric(adj_pval < my_alpha))
  
  batchbh_res <- data %>%
    select(id, batch, pval) %>%
    BatchBH(gammai = c(0.5, 0.5, rep(0, n_batches - 2))) #aggressive investing schedule
  
  batchstbh_res <- data %>%
    select(id, batch, pval) %>%
    BatchStBH(gammai = c(0.5, 0.5, rep(0, n_batches - 2))) #aggressive investing schedule
  
  batchprds_res <- data %>%
    select(id, batch, pval) %>%
    BatchPRDS() #default investing schedule
  
  out <- data %>%
    left_join(batchbh_res %>% select(id, batch, R), by = c("id", "batch")) %>%
    rename(R_on_bh = R) %>%
    left_join(batchstbh_res %>% select(id, batch, R), by = c("id", "batch")) %>%
    rename(R_on_stbh = R) %>%
    left_join(batchprds_res %>% select(id, batch, R), by = c("id", "batch")) %>%
    rename(R_on_prds = R) %>%
    left_join(offline_bh_res %>% select(id, batch, R_off_bh), 
              by = c("id", "batch")) %>%
    left_join(offline_stbh_res %>% select(id, batch, R_off_stbh), 
              by = c("id", "batch")) %>%
    left_join(pooled_offline_bh_res %>% select(id, batch, R_pool_bh),
              by = c("id", "batch")) %>%
    left_join(pooled_offline_stbh_res %>% select(id, batch, R_pool_stbh),
              by = c("id", "batch")) %>%
    mutate(R_0.025 = as.numeric(pval < 0.025),
           R_0.05 = as.numeric(pval < 0.05),
           R_0.1 = as.numeric(pval < 0.1))
  out
}

all_fdr_res <- list()
for (i in seq_along(1:length(all_rnaseq_de_res))) {
  all_fdr_res[[i]] <- apply_all_fdr_procedures(all_rnaseq_de_res[[i]])
}

# Calculate metrics (power + fdr) ----
calculate_metrics <- function(data, .pi) {
  
  #calculate metrics for each data element (defined by pi)
  R_0.025 <- data$R_0.025
  R_0.05 <- data$R_0.05
  R_0.1 <- data$R_0.1
  R_off_bh <- data$R_off_bh
  R_off_stbh <- data$R_off_stbh
  R_pool_bh <- data$R_pool_bh
  R_pool_stbh <- data$R_pool_stbh
  R_on_bh <- data$R_on_bh
  R_on_stbh <- data$R_on_stbh
  R_on_prds <- data$R_on_prds
  truth <- data$truth
  
  positives <- sum(truth, na.rm = TRUE)
  negatives <- sum(truth == 0, na.rm = TRUE)
  tp_0.025 <- sum(R_0.025 & truth, na.rm = TRUE)
  fp_0.025 <- sum(R_0.025 & truth == 0, na.rm = TRUE)
  tn_0.025 <- sum(R_0.025 == 0 & truth == 0, na.rm = TRUE)
  fn_0.025 <- sum(R_0.025 == 0 & truth, na.rm = TRUE)
  power_0.025 <- tp_0.025/positives
  fdr_0.025 <- fp_0.025/sum(R_0.025, na.rm = TRUE)
  tpr_0.025 <- tp_0.025/(tp_0.025 + fn_0.025)
  tnr_0.025 <- tn_0.025/(tn_0.025 + fp_0.025)
  ppv_0.025 <- tp_0.025/(tp_0.025 + fp_0.025)
  
  tp_0.05 <- sum(R_0.05 & truth, na.rm = TRUE)
  fp_0.05 <- sum(R_0.05 & truth == 0, na.rm = TRUE)
  tn_0.05 <- sum(R_0.05 == 0 & truth == 0, na.rm = TRUE)
  fn_0.05 <- sum(R_0.05 == 0 & truth, na.rm = TRUE)
  power_0.05 <- tp_0.05/positives
  fdr_0.05 <- fp_0.05/sum(R_0.05, na.rm = TRUE)
  tpr_0.05 <- tp_0.05/(tp_0.05 + fn_0.05)
  tnr_0.05 <- tn_0.05/(tn_0.05 + fp_0.05)
  ppv_0.05 <- tp_0.05/(tp_0.05 + fp_0.05)
  
  tp_0.1 <- sum(R_0.1 & truth, na.rm = TRUE)
  fp_0.1 <- sum(R_0.1 & truth == 0, na.rm = TRUE)
  tn_0.1 <- sum(R_0.1 == 0 & truth == 0, na.rm = TRUE)
  fn_0.1 <- sum(R_0.1 == 0 & truth, na.rm = TRUE)
  power_0.1 <- tp_0.1/positives
  fdr_0.1 <- fp_0.1/sum(R_0.1, na.rm = TRUE)
  tpr_0.1 <- tp_0.1/(tp_0.1 + fn_0.1)
  tnr_0.1 <- tn_0.1/(tn_0.1 + fp_0.1)
  ppv_0.1 <- tp_0.1/(tp_0.1 + fp_0.1)
  
  tp_off_bh <- sum(R_off_bh & truth, na.rm = TRUE)
  fp_off_bh <- sum(R_off_bh & truth == 0, na.rm = TRUE)
  tn_off_bh <- sum(R_off_bh == 0 & truth == 0, na.rm = TRUE)
  fn_off_bh <- sum(R_off_bh == 0 & truth, na.rm = TRUE)
  power_off_bh <- tp_off_bh/positives
  fdr_off_bh <- fp_off_bh/sum(R_off_bh, na.rm = TRUE)
  tpr_off_bh <- tp_off_bh/(tp_off_bh + fn_off_bh)
  tnr_off_bh <- tn_off_bh/(tn_off_bh + fp_off_bh)
  ppv_off_bh <- tp_off_bh/(tp_off_bh + fp_off_bh)
  
  tp_off_stbh <- sum(R_off_stbh & truth, na.rm = TRUE)
  fp_off_stbh <- sum(R_off_stbh & truth == 0, na.rm = TRUE)
  tn_off_stbh <- sum(R_off_stbh == 0 & truth == 0, na.rm = TRUE)
  fn_off_stbh <- sum(R_off_stbh == 0 & truth, na.rm = TRUE)
  power_off_stbh <- tp_off_stbh/positives
  fdr_off_stbh <- fp_off_stbh/sum(R_off_stbh, na.rm = TRUE)
  tpr_off_stbh <- tp_off_stbh/(tp_off_stbh + fn_off_stbh)
  tnr_off_stbh <- tn_off_stbh/(tn_off_stbh + fp_off_stbh)
  ppv_off_stbh <- tp_off_stbh/(tp_off_stbh + fp_off_stbh)
  
  tp_pool_bh <- sum(R_pool_bh & truth, na.rm = TRUE)
  fp_pool_bh <- sum(R_pool_bh & truth == 0, na.rm = TRUE)
  tn_pool_bh <- sum(R_pool_bh == 0 & truth == 0, na.rm = TRUE)
  fn_pool_bh <- sum(R_pool_bh == 0 & truth, na.rm = TRUE)
  power_pool_bh <- tp_pool_bh/positives
  fdr_pool_bh <- fp_pool_bh/sum(R_pool_bh, na.rm = TRUE)
  tpr_pool_bh <- tp_pool_bh/(tp_pool_bh + fn_pool_bh)
  tnr_pool_bh <- tn_pool_bh/(tn_pool_bh + fp_pool_bh)
  ppv_pool_bh <- tp_pool_bh/(tp_pool_bh + fp_pool_bh)
  
  tp_pool_stbh <- sum(R_pool_stbh & truth, na.rm = TRUE)
  fp_pool_stbh <- sum(R_pool_stbh & truth == 0, na.rm = TRUE)
  tn_pool_stbh <- sum(R_pool_stbh == 0 & truth == 0, na.rm = TRUE)
  fn_pool_stbh <- sum(R_pool_stbh == 0 & truth, na.rm = TRUE)
  power_pool_stbh <- tp_pool_stbh/positives
  fdr_pool_stbh <- fp_pool_stbh/sum(R_pool_stbh, na.rm = TRUE)
  tpr_pool_stbh <- tp_pool_stbh/(tp_pool_stbh + fn_pool_stbh)
  tnr_pool_stbh <- tn_pool_stbh/(tn_pool_stbh + fp_pool_stbh)
  ppv_pool_stbh <- tp_pool_stbh/(tp_pool_stbh + fp_pool_stbh)
  
  tp_on_bh <- sum(R_on_bh & truth, na.rm = TRUE)
  fp_on_bh <- sum(R_on_bh & truth == 0, na.rm = TRUE)
  tn_on_bh <- sum(R_on_bh == 0 & truth == 0, na.rm = TRUE)
  fn_on_bh <- sum(R_on_bh == 0 & truth, na.rm = TRUE)
  power_on_bh <- tp_on_bh/positives
  fdr_on_bh <- fp_on_bh/sum(R_on_bh, na.rm = TRUE)
  tpr_on_bh <- tp_on_bh/(tp_on_bh + fn_on_bh)
  tnr_on_bh <- tn_on_bh/(tn_on_bh + fp_on_bh)
  ppv_on_bh <- tp_on_bh/(tp_on_bh + fp_on_bh)
  
  tp_on_stbh <- sum(R_on_stbh & truth, na.rm = TRUE)
  fp_on_stbh <- sum(R_on_stbh & truth == 0, na.rm = TRUE)
  tn_on_stbh <- sum(R_on_stbh == 0 & truth == 0, na.rm = TRUE)
  fn_on_stbh <- sum(R_on_stbh == 0 & truth, na.rm = TRUE)
  power_on_stbh <- tp_on_stbh/positives
  fdr_on_stbh <- fp_on_stbh/sum(R_on_stbh, na.rm = TRUE)
  tpr_on_stbh <- tp_on_stbh/(tp_on_stbh + fn_on_stbh)
  tnr_on_stbh <- tn_on_stbh/(tn_on_stbh + fp_on_stbh)
  ppv_on_stbh <- tp_on_stbh/(tp_on_stbh + fp_on_stbh)
  
  tp_on_prds <- sum(R_on_prds & truth, na.rm = TRUE)
  fp_on_prds <- sum(R_on_prds & truth == 0, na.rm = TRUE)
  tn_on_prds <- sum(R_on_prds == 0 & truth == 0, na.rm = TRUE)
  fn_on_prds <- sum(R_on_prds == 0 & truth, na.rm = TRUE)
  power_on_prds <- tp_on_prds/positives
  fdr_on_prds <- fp_on_prds/sum(R_on_prds, na.rm = TRUE)
  tpr_on_prds <- tp_on_prds/(tp_on_prds + fn_on_prds)
  tnr_on_prds <- tn_on_prds/(tn_on_prds + fp_on_prds)
  ppv_on_prds <- tp_on_prds/(tp_on_prds + fp_on_prds)
  
  tibble(
    pi = .pi,
    tp_0.025 = tp_0.025,
    fp_0.025 = fp_0.025,
    tn_0.025 = tn_0.025,
    fn_0.025 = fn_0.025,
    power_0.025 = power_0.025,
    fdr_0.025 = fdr_0.025,
    tpr_0.025 = tpr_0.025,
    tnr_0.025 = tnr_0.025,
    ppv_0.025 = ppv_0.025,
    tp_0.05 = tp_0.05,
    fp_0.05 = fp_0.05,
    tn_0.05 = tn_0.05,
    fn_0.05 = fn_0.05,
    power_0.05 = power_0.05,
    fdr_0.05 = fdr_0.05,
    tpr_0.05 = tpr_0.05,
    tnr_0.05 = tnr_0.05,
    ppv_0.05 = ppv_0.05,
    tp_0.1 = tp_0.1,
    fp_0.1 = fp_0.1,
    tn_0.1 = tn_0.1,
    fn_0.1 = fn_0.1,
    power_0.1 = power_0.1,
    fdr_0.1 = fdr_0.1,
    tpr_0.1 = tpr_0.1,
    tnr_0.1 = tnr_0.1,
    ppv_0.1 = ppv_0.1,
    tp_off_bh = tp_off_bh,
    fp_off_bh = fp_off_bh,
    tn_off_bh = tn_off_bh,
    fn_off_bh = fn_off_bh,
    power_off_bh = power_off_bh,
    fdr_off_bh = fdr_off_bh,
    tpr_off_bh = tpr_off_bh,
    tnr_off_bh = tnr_off_bh,
    ppv_off_bh = ppv_off_bh,
    tp_off_stbh = tp_off_stbh,
    fp_off_stbh = fp_off_stbh,
    tn_off_stbh = tn_off_stbh,
    fn_off_stbh = fn_off_stbh,
    power_off_stbh = power_off_stbh,
    fdr_off_stbh = fdr_off_stbh,
    tpr_off_stbh = tpr_off_stbh,
    tnr_off_stbh = tnr_off_stbh,
    ppv_off_stbh = ppv_off_stbh,
    tp_pool_bh = tp_pool_bh,
    fp_pool_bh = fp_pool_bh,
    tn_pool_bh = tn_pool_bh,
    fn_pool_bh = fn_pool_bh,
    power_pool_bh = power_pool_bh,
    fdr_pool_bh = fdr_pool_bh,
    tpr_pool_bh = tpr_pool_bh,
    tnr_pool_bh = tnr_pool_bh,
    ppv_pool_bh = ppv_pool_bh,
    tp_pool_stbh = tp_pool_stbh,
    fp_pool_stbh = fp_pool_stbh,
    tn_pool_stbh = tn_pool_stbh,
    fn_pool_stbh = fn_pool_stbh,
    power_pool_stbh = power_pool_stbh,
    fdr_pool_stbh = fdr_pool_stbh,
    tpr_pool_stbh = tpr_pool_stbh,
    tnr_pool_stbh = tnr_pool_stbh,
    ppv_pool_stbh = ppv_pool_stbh,
    tp_on_bh = tp_on_bh,
    fp_on_bh = fp_on_bh,
    tn_on_bh = tn_on_bh,
    fn_on_bh = fn_on_bh,
    power_on_bh = power_on_bh,
    fdr_on_bh = fdr_on_bh,
    tpr_on_bh = tpr_on_bh,
    tnr_on_bh = tnr_on_bh,
    ppv_on_bh = ppv_on_bh,
    tp_on_stbh = tp_on_stbh,
    fp_on_stbh = fp_on_stbh,
    tn_on_stbh = tn_on_stbh,
    fn_on_stbh = fn_on_stbh,
    power_on_stbh = power_on_stbh,
    fdr_on_stbh = fdr_on_stbh,
    tpr_on_stbh = tpr_on_stbh,
    tnr_on_stbh = tnr_on_stbh,
    ppv_on_stbh = ppv_on_stbh,
    tp_on_prds = tp_on_prds,
    fp_on_prds = fp_on_prds,
    tn_on_prds = tn_on_prds,
    fn_on_prds = fn_on_prds,
    power_on_prds = power_on_prds,
    fdr_on_prds = fdr_on_prds,
    tpr_on_prds = tpr_on_prds,
    tnr_on_prds = tnr_on_prds,
    ppv_on_prds = ppv_on_prds
  )
}

metrics_out <- list()

for (i in seq_along(1:length(all_fdr_res))) {
  metrics_out[[i]] <- calculate_metrics(all_fdr_res[[i]], as.numeric(names(pi_list)[i]))
}

all_metrics_out <- bind_rows(metrics_out)

saveRDS(all_metrics_out,
        paste0(output_dir, "change_pi_run_reorder1.5_", task_id, ".rds"))