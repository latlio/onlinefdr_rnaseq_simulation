# Functions to simulate RNASeq dataset

# ARCHIVE ----
simulate_simple_rnaseq <- function(n_genes,
                                   raw_count_exp_rate = 0.01,
                                   prop_DE = 0.05,
                                   logfc_sd = 2,
                                   sample_nb_means,
                                   n_samples_per_group,
                                   dispersion_parameter = 0.2) {
  raw_counts <- rexp(n_genes, rate = raw_count_exp_rate)
  is_de <- runif(n_genes) < prop_DE
  lfc <- rnorm(n_genes, sd = logfc_sd)
  de_counts_up <- ifelse(is_de, raw_counts * 2^(lfc/2), raw_counts)
  de_counts_down <- ifelse(is_de, raw_counts * 2^(-lfc/2), raw_counts)
  
  total_samples <- c(rep("control", n_samples_per_group),
                     rep("exp", n_samples_per_group))
  if(length(sample_nb_means) != length(total_samples)) {
    stop("Please provide a vector of estimated means whose length matches that 
         of the total number of samples.")
  }
  
  count_matrix <- t(sapply(seq_len(n_genes), function(i)
    sapply(seq(length(total_samples)), function(j)
      rnbinom(1, mu = sample_nb_means[j] * ifelse(total_samples[j]=="control", 
                                                  de_counts_up[i],
                                                  de_counts_down[i]),
               size =  dispersion_parameter))))
  colnames(count_matrix) <- paste0(total_samples, rep(seq(n_samples_per_group), 
                                                      2))
  rownames(count_matrix) <- paste("gene", seq_len(n_genes),
                                  ifelse(is_de, "T", "F" ), sep="_" )
  count_matrix
}

swap_indices <- function(x,i,j) {x[c(i,j)] <- x[c(j,i)]; x} 

create_per_gene_effect_sizes <- function(n_genes,
                                         prop_diff_exp,
                                         desired_effect_size,
                                         noise_factor = 5) {
  #' n_genes: number of genes
  #' prop_diff_exp: numeric, insert a proportion (p), which will be used
  #' to assign the first p*n_genes genes as differentially expressed
  #' desired_effect_size: numeric, insert a value greater than 0
  #' noise_factor: numeric, number of genes you'd like to swap from being differentially
  #' expressed to non-differentially expressed and vice versa, calculated as
  #' n_genes*prop_diff_exp/noise_factor
  #' Output: vector of length n_genes that specifies the effect_sizes for each gene,
  #' with the first p*n_genes being diff. exp. at an effect size specified by 
  #' desired_effect_size and the remaining n_genes - p*n_genes genes with an
  #' effect size of 1
  #' Adding noise factor will randomly select some of the genes that are supposed to be 
  #' diff. exp. and convert them to null effect size, which aims to simulate natural
  #' variability of certain DE genes not being picked up as DE (e.g. low signal)
  
  out <- c(rep(desired_effect_size, each = n_genes*prop_diff_exp), 
           rep(1, each = n_genes - (n_genes*prop_diff_exp)))
  
  if(is.null(noise_factor)) {
    return(out)
  } else {
    #' swap certain genes to be differentially expressed to simulate natural
    #' variability
    if(n_genes*prop_diff_exp/noise_factor >= n_genes*prop_diff_exp){
      stop("The number of genes you swap must be
  less than the number of differentially expressed genes.")
    }
    swapped_out <- swap_indices(out, 
                                sample.int(n_genes*prop_diff_exp, n_genes*prop_diff_exp/noise_factor),
                                sample.int(n_genes, n_genes*prop_diff_exp/noise_factor))
    return(swapped_out)
  }
}

# CURRENT ----
create_per_gene_effect_sizes2 <- function(n_genes, pi, fc) {
  c(rep(2^fc, n_genes*pi/2), rep(1/(2^fc), n_genes*pi/2), rep(1, n_genes-(n_genes*pi)))
}

#differs from tidyde::run_sample_limma_de in that it doesn't require metadata input
run_limma_de <- function(n_samples, 
                         sim_counts) {
  #this function is meant to be run in a map
  # assumes same samples between treatment groups
  sim_meta <- tibble(sample = colnames(sim_counts),
                     condition = c(rep("control", n_samples), 
                                   rep("experimental", n_samples)))
  sim_design <- check_sample_names(sim_counts,
                                   NULL,
                                   sim_meta, 
                                   sample) %>%
    purrr::pluck("meta") %>%
    make_design_matrix(., c("condition"))
  
  sim_id <- paste0("g", seq(nrow(sim_counts)))
  
  sim_de_res <- check_sample_names(sim_counts, 
                                   NULL, 
                                   sim_meta, 
                                   sample) %>%
    purrr::pluck("mod_count") %>%
    filter_genes(., sim_id, "edgeR") %>%
    make_voom(., sim_design) %>%
    model_limma() %>%
    make_contrasts(design_matrix = sim_design, 
                   "conditioncontrol",
                   "conditionexperimental") %>%
    model_bayes()
  
  tidy.marray.lm(sim_de_res)
}

apply_fdr <- function(data, n_batch, offline_method = "BH", alpha = 0.05) {
  
  pooled_data <- data %>%
    group_by(batch) %>%
    mutate(uncorrected_0.1 = as.numeric(pval < 0.1),
           uncorrected_0.05 = as.numeric(pval < 0.05),
           uncorrected_0.025 = as.numeric(pval < 0.025),
           offline_naive_adj_pval = p.adjust(pval, method = offline_method),
           offline_naive_R = as.numeric(offline_naive_adj_pval < alpha)) %>%
    ungroup() 
  
  #this creates incrementally growing families at each "time point"
  pooled_data_families <- list()
  for(i in seq_along(1:n_batch)) {
    pooled_data_families[[i]] <- pooled_data %>%
      filter(batch %in% c(1:i))
  }
  
  #apply offline pooled procedures to each pooled data family
  pooled_data_families <- map(pooled_data_families, ~.x %>% 
                                mutate(offline_pooled_adj_pval = p.adjust(pval, method = offline_method),
                                       offline_pooled_R = as.numeric(offline_pooled_adj_pval < alpha)))
  
  batchbh_data <- onlineFDR::BatchBH(pooled_data %>%
                                       select(id, batch, pval),
                                     alpha = alpha,
                                     gammai = c(0.5, 0.5, rep(0, n_batch - 2))) %>%
    rename(batchbh_R = R) %>%
    select(id, batch, batchbh_R)
  
  batchprds_data <- onlineFDR::BatchPRDS(pooled_data %>%
                                           select(id, batch, pval),
                                         alpha = alpha) %>%
    rename(batchprds_R = R) %>%
    select(id, batch, batchprds_R)
  
  batchstbh_data <- onlineFDR::BatchStBH(pooled_data %>%
                                           select(id, batch, pval),
                                         alpha = alpha,
                                         gammai = c(0.5, 0.5, rep(0, n_batch - 2))) %>%
    rename(batchstbh_R = R) %>%
    select(id, batch, batchstbh_R)
  
  #join online results and truth to each pooled_data_family for subsequent analysis
  pooled_data_families <- map(pooled_data_families, ~.x %>%
                                left_join(batchbh_data, by = c("id", "batch")) %>%
                                left_join(batchprds_data, by = c("id", "batch")) %>%
                                left_join(batchstbh_data, by = c("id", "batch")))
  
  pooled_data_families
}

calculate_incremental_power <- function(df) {
  df %>%
    group_by(batch) %>%
    summarise(numerator_uncorrected_0.1 = sum(uncorrected_0.1 & truth),
              numerator_uncorrected_0.05 = sum(uncorrected_0.05 & truth),
              numerator_uncorrected_0.025 = sum(uncorrected_0.025 & truth),
              numerator_offline_naive = sum(offline_naive_R & truth),
              numerator_batchbh = sum(batchbh_R & truth),
              numerator_batchprds = sum(batchprds_R & truth),
              numerator_batchstbh = sum(batchstbh_R & truth),
              denominator = sum(truth)) %>%
    transmute(batch, 
              power_uncorrected_0.1 = numerator_uncorrected_0.1/denominator,
              power_uncorrected_0.05 = numerator_uncorrected_0.05/denominator,
              power_uncorrected_0.025 = numerator_uncorrected_0.025/denominator,
              power_offline_naive = numerator_offline_naive/denominator,
              power_batchbh = cumsum(numerator_batchbh)/cumsum(denominator),
              power_batchprds = cumsum(numerator_batchprds)/cumsum(denominator),
              power_batchstbh = cumsum(numerator_batchstbh)/cumsum(denominator)) %>%
    ungroup()
}

calculate_offline_pooled_power <- function(df) {
  df %>%
    summarize(power_offline_pooled = sum(offline_pooled_R & truth) / sum(truth))
}

calculate_incremental_fdp <- function(df) {
  df %>%
    group_by(batch) %>%
    summarise(numerator_uncorrected_0.1 = sum(uncorrected_0.1 & truth == 0),
              numerator_uncorrected_0.05 = sum(uncorrected_0.05 & truth == 0),
              numerator_uncorrected_0.025 = sum(uncorrected_0.025 & truth == 0),
              numerator_offline_naive = sum(offline_naive_R & truth == 0),
              numerator_offline_pooled = sum(offline_pooled_R & truth == 0), 
              numerator_batchbh = sum(batchbh_R & truth == 0),
              numerator_batchprds = sum(batchprds_R & truth == 0),
              numerator_batchstbh = sum(batchstbh_R & truth == 0),
              denominator_uncorrected_0.1 = sum(uncorrected_0.1),
              denominator_uncorrected_0.05 = sum(uncorrected_0.05),
              denominator_uncorrected_0.025 = sum(uncorrected_0.025),
              denominator_offline_naive = sum(offline_naive_R),
              denominator_offline_pooled = sum(offline_pooled_R),
              denominator_batchbh = sum(batchbh_R),
              denominator_batchprds = sum(batchprds_R),
              denominator_batchstbh = sum(batchstbh_R)) %>%
    transmute(batch, 
              fdp_uncorrected_0.1 = numerator_uncorrected_0.1/denominator_uncorrected_0.1,
              fdp_uncorrected_0.05 = numerator_uncorrected_0.05/denominator_uncorrected_0.05,
              fdp_uncorrected_0.025 = numerator_uncorrected_0.025/denominator_uncorrected_0.025,
              fdp_offline_naive = numerator_offline_naive/denominator_offline_naive,
              fdp_offline_pooled = cumsum(numerator_offline_pooled)/cumsum(denominator_offline_pooled),
              fdp_batchbh = cumsum(numerator_batchbh)/cumsum(denominator_batchbh),
              fdp_batchprds = cumsum(numerator_batchprds)/cumsum(denominator_batchprds),
              fdp_batchstbh = cumsum(numerator_batchstbh)/cumsum(denominator_batchstbh)) %>%
    ungroup()
}

calculate_incremental_discrimination <- function(df) {
  df %>%
    group_by(batch) %>%
    summarize(TP_offline_pooled = sum(offline_pooled_R & truth),
              TN_offline_pooled = sum(offline_pooled_R == 0 & truth == 0),
              FP_offline_pooled = sum(offline_pooled_R & truth == 0),
              FN_offline_pooled = sum(offline_pooled_R == 0 & truth),
              # TP_offline_naive = sum(offline_naive_R & truth),
              # TN_offline_naive = sum(offline_naive_R == 0 & truth == 0),
              # FP_offline_naive = sum(offline_naive_R & truth == 0),
              # FN_offline_naive = sum(offline_naive_R == 0 & truth),
              TP_batchbh = sum(batchbh_R & truth),
              TN_batchbh = sum(batchbh_R == 0 & truth == 0),
              FP_batchbh = sum(batchbh_R & truth == 0),
              FN_batchbh = sum(batchbh_R == 0 & truth),
              TP_batchprds = sum(batchprds_R & truth),
              TN_batchprds = sum(batchprds_R == 0 & truth == 0),
              FP_batchprds = sum(batchprds_R & truth == 0),
              FN_batchprds = sum(batchprds_R == 0 & truth),
              TP_batchstbh = sum(batchstbh_R & truth),
              TN_batchstbh = sum(batchstbh_R == 0 & truth == 0),
              FP_batchstbh = sum(batchstbh_R & truth == 0),
              FN_batchstbh = sum(batchstbh_R == 0 & truth)) %>%
    transmute(batch,
              TP_offline_pooled_cum = cumsum(TP_offline_pooled),
              TN_offline_pooled_cum = cumsum(TN_offline_pooled),
              FP_offline_pooled_cum = cumsum(FP_offline_pooled),
              FN_offline_pooled_cum = cumsum(FN_offline_pooled),
              # TP_offline_naive,
              # TN_offline_naive,
              # FP_offline_naive,
              # FN_offline_naive,
              TP_batchbh_cum = cumsum(TP_batchbh),
              TN_batchbh_cum = cumsum(TN_batchbh),
              FP_batchbh_cum = cumsum(FP_batchbh),
              FN_batchbh_cum = cumsum(FN_batchbh),
              TP_batchprds_cum = cumsum(TP_batchprds),
              TN_batchprds_cum = cumsum(TN_batchprds),
              FP_batchprds_cum = cumsum(FP_batchprds),
              FN_batchprds_cum = cumsum(FN_batchprds),
              TP_batchstbh_cum = cumsum(TP_batchstbh),
              TN_batchstbh_cum = cumsum(TN_batchstbh),
              FP_batchstbh_cum = cumsum(FP_batchstbh),
              FN_batchstbh_cum = cumsum(FN_batchstbh)) %>%
    ungroup()
}

#' Functions for simulating RNAseq DE data
#' 
#' Check sample names
#'
#' `check_sample_names()` is a simple quality control step that verifies whether
#' the column names in the count matrix match with a user-defined metadata column.
#' In order for a match to occur, the value levels of the column names and those of the
#' user-defined metadata column need to be identical
#' (e.g. `setdiff(colnames(my_count_matrix), metadata$my_column)` should be 0),
#' and the order in which the values appear need to be identical
#' (e.g. `identical(colnames(my_count_matrix), metadata$my_column)` should be TRUE).
#'
#' @param count_df cleaned dataframe of counts, rows should be gene IDs,
#' columns should be samples, cells should only contain counts
#' @param cols_to_remove vector of column numbers that do not correspond to a sample,
#' necessary to identify for downstream functions
#' @param metadata cleaned metadata for RNAseq data
#' @param metadata_var column of sample identifier that user expects to match with
#' count_matrix
#'
#' @return a `list` with the following components:
#' \item{old_count}{the original count dataframe supplied}
#' \item{mod_count}{the pure count dataframe (no other columns)}
#' \item{meta}{sorted metadata (if necessary), otherwise the supplied
#' metadata is returned with console message output of the quality control check.}
#'
#' @details The proportion of zeros in the original count data is also printed to the
#' console. For count data that has a medium to high proportion of zeros,
#' \code{\link[edgeR]{voomLmFit}} is recommended. Otherwise
#' \code{\link[edgeR]{voom}} followed by \code{\link[edgeR]{lmFit}} is
#' recommended.
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#' mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' check_sample_names(counts, c(1,2), meta, FileName)
#'
check_sample_names <- function(count_df,
                               cols_to_remove,
                               metadata,
                               metadata_var) {
  # coerce to dataframe
  if(!is.data.frame(count_df)) {
    message("Your input data was coerced into a data frame.")
    count_df <- as.data.frame(count_df)
  }
  
  #output count dataframe with purely counts
  if(is.null(cols_to_remove)) {
    count_df_mod <- count_df
  } else {
    count_df_mod <- count_df[,-cols_to_remove]
    if (
      identical(
        setdiff(colnames(count_df)[-cols_to_remove],
                metadata %>% dplyr::select({{metadata_var}}) %>% dplyr::pull()),
        character(0)
      )
    ){
      message("The column names of the count matrix and the unique sample ID values are correctly specified.")
      
      if(
        identical(
          colnames(count_df)[-cols_to_remove],
          metadata %>% dplyr::select({{metadata_var}}) %>% dplyr::pull()
        )
      ) {
        message("The order is also correct. You can safely proceed with the remaining analysis steps.")
      } else {
        message("The data was ordered incorrectly. The metadata has been reordered.")
        metadata <- metadata[
          match(
            colnames(count_df)[-cols_to_remove],
            metadata %>% dplyr::select({{metadata_var}}) %>% dplyr::pull()
          ),]
      }
    } else {
      stop("Please specify a different variable or check that the values in the metadata are written correctly.")
    }
  }
  
  cat("The proportion of zeroes in your count data is ",
      sum(count_df == 0)/(ncol(count_df) * nrow(count_df)))
  
  list(
    old_count = count_df,
    mod_count = count_df_mod,
    meta = metadata
  )
}

#' Create a design matrix
#'
#' `make_design_matrix()` creates a model matrix for your DE analysis.
#' See \code{\link[stats]{model.matrix}} for further details on what a design
#' (or model) matrix is.
#'
#' @param metadata cleaned metadata for RNAseq data
#' @param vars a character vector of variables to include in the model
#'
#' @details The order in which you specify your variables will affect which
#' variable is dummy coded as the reference variable.
#'
#' @return a `tbl` of the design matrix
#'
#' @export
#'
#' @examples
#' make_design_matrix(metadata, "CellType")
#'
#' make_design_matrix(metadata, c("CellType", "Status"))
#'
#' # Results in a different design matrix
#' make_design_matrix(metadata, c("Status", "CellType"))
#'
#' # In a pipeline
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' my_design <- check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("meta") %>%
#'   make_design_matrix(., c("CellType"))

make_design_matrix <- function(metadata, vars) {
  formula <- as.formula(paste0("~ 0 + ", paste(vars, collapse = "+")))
  modelr::model_matrix(metadata,
                       formula)
}

#'
#' Filter lowly expressed genes
#'
#' `filter_genes()` is a wrapper function for several filtering methods.
#'
#' @param count_df preprocessed dataframe of pure counts
#' @param id vector of gene IDs
#' @param filter_method Either `edgeR`, `samplenr`, or `cpm`
#' @param min_samples minimum number of samples
#' @param min_cpm minimum cpm
#' @param ... additional arguments to `filterByExpr()`
#'
#' @details I encourage users to exercise caution before using this filter function.
#' Oftentimes, the filtering step should be specific to the sequencing experiment.
#' The `edgeR` option is a wrapper for `edgeR::filterByExpr()`.
#' The `samplenr` option filters out genes across sample whose counts are lower
#' 2*number_of_samples
#' The `cpm` option filters out genes whose rowsums (excluding cells lower
#' than `min_cpm`) are less than number_of_samples/min_samples
#'
#' @return a `list` (`DGEList`) with the following components:
#' \item{counts}{a vector of the filtered counts}
#' \item{samples}{a dataframe containing the library sizes and the normalization
#' factors}
#' \item{genes}{a dataframe containing the gene IDs}
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' # this step may differ depending on how your data is formatted
#' id <- as.character(counts$EntrezGeneID)
#'
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR")

filter_genes <- function(count_df,
                         id,
                         filter_method,
                         min_samples = 10,
                         min_cpm = 0.25,
                         ...) {
  
  filter_by_edgeR <- function(dge, id) {
    counts <- dge$counts
    #rownames is not "tidy" but it may be easier to work with
    rownames(counts) <- id
    genes_to_keep <- edgeR::filterByExpr(dge, ...)
    counts <- counts[genes_to_keep,]
    genes <- rownames(counts)
    dge_filtered <- edgeR::DGEList(counts = counts,
                                   genes = genes)
    dge_filtered
  }
  
  filter_by_samplenr <- function(dge, id) {
    counts <- dge$counts
    rownames(counts) <- id
    genes_to_keep <- rowSums(counts) >= 2*ncol(counts)
    counts <- counts[genes_to_keep,]
    genes <- rownames(counts)
    dge_filtered <- edgeR::DGEList(counts = counts,
                                   genes = genes)
    dge_filtered
  }
  
  filter_by_cpm <- function(dge,
                            id,
                            min_samples,
                            min_cpm) {
    counts <- dge$counts
    rownames(counts) <- id
    tmp <- edgeR::cpm(counts) %>% as.data.frame()
    
    genes_to_keep <- tmp %>%
      tibble::rownames_to_column("id") %>%
      rowwise() %>%
      mutate(cond = {ind <- c_across(2:last_col())
      sum(ind >= min_cpm) >= ncol(tmp)/min_samples}) %>%
      ungroup() %>%
      tibble::column_to_rownames("id") %>%
      filter(cond)
    
    counts <- counts[rownames(genes_to_keep),]
    genes <- rownames(counts)
    dge_filtered <- edgeR::DGEList(counts = counts,
                                   genes = genes)
    dge_filtered
  }
  
  dge <- edgeR::DGEList(count_df)
  dge_filtered <- switch(filter_method,
                         edgeR = filter_by_edgeR(dge, id),
                         samplenr = filter_by_samplenr(dge, id),
                         cpm = filter_by_cpm(dge,
                                             id,
                                             min_samples,
                                             min_cpm))
  dgelist_w_normfactors <- calcNormFactors(dge_filtered)
  dgelist_w_normfactors
}

#' Convert count data to voom
#'
#' `make_voom()` is a wrapper function for `voom()`
#'
#' @param .dge a `list` or a `DGElist` object
#' @param design_matrix a design matrix with rows corresponding to samples
#' and columns to coefficients to be estimated
#' @param .f limma::voom
#' @param ... additional arguments passed to `.f`
#'
#' @details Please refer to \code{\link[limma]{voom}} for more information.
#'
#' @return a `list` with the following components:
#' \item{E}{a numeric matrix of normalized expression values on the log2 scale}
#' \item{weights}{numeric matrix of inverse variance weights}
#' \item{design}{design matrix}
#' \item{lib.size}{numeric vector of total normalized library sizes}
#' \item{genes}{dataframe of gene annotation extracted from `counts`}
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' id <- as.character(counts$EntrezGeneID)
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %>%
#'   make_voom(., my_design)

make_voom <- function(.dge, design_matrix, .f = limma::voom, ...) {
  .args <- rlang::enexprs(...)
  
  rlang::eval_tidy(rlang::expr(.f(counts = .dge,
                                  design = design_matrix,
                                  !!! .args)))
}

#' Model limma
#'
#' `model_limma()` is a wrapper function for `lmFit`. Fits a linear model for
#' each gene.
#'
#' @param voom a voom object
#' @param .f limma::lmFit
#' @param ... additional arguments to .f
#'
#' @details Please refer to \code{\link[limma]{lmFit}} for more information.
#'
#' @return a `list` (`MArrayLM`) object containing the results of the fit
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' id <- as.character(counts$EntrezGeneID)
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %>%
#'   make_voom(., my_design) %>%
#'   model_limma()

model_limma <- function(.data, .f = limma::lmFit, ...) {
  .args <- rlang::enexprs(...)
  
  rlang::eval_tidy(rlang::expr(.f(object = .data,
                                  !!! .args)))
}

#' Make contrasts
#'
#' `make_contrasts()` computes estimated coefficients and standard errors
#' for a given set of contrasts
#'
#' @param .fit an MArrayLM object or list object. Must contain components
#' `coefficients` and `stdev.unscaled`
#' @param design_matrix a design matrix with rows corresponding to samples
#' and columns to coefficients to be estimated
#' @param ... the names of the variables which you like to compare
#'
#' @details Please refer to \code{\link[limma]{contrasts.fit}} for more information.
#'
#' @return a `list` object of the same class as .fit
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' id <- as.character(counts$EntrezGeneID)
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %>%
#'   make_voom(., my_design) %>%
#'   model_limma() %>%
#'   make_contrasts(Statuspregnant, Statusvirgin)
make_contrasts <- function(.fit, design_matrix, ...) {
  contrast_components <- rlang::enexprs(...)
  my_contrast <- paste(contrast_components, collapse = " - ")
  contrast_matrix <- limma::makeContrasts(
    contrasts = c(my_contrast),
    levels = colnames(design_matrix)
  )
  contrasts.fit(.fit, contrast_matrix)
}

#' Model differential expression
#'
#' `model_bayes()` performs an empirical Bayes fit
#'
#' @param .fit an MArrayLM object produced by `model_limma()` or
#' `make_contrasts()`
#' @param .f limma::eBayes
#' @param ... additional arguments to .f
#'
#' @details Please refer to \code{\link[limma]{eBayes}} for more information.
#'
#' @return a `list` (`MArrayLM`) object containing the results of the fit
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' id <- as.character(counts$EntrezGeneID)
#' check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %>%
#'   make_voom(., my_design) %>%
#'   model_limma() %>%
#'   make_contrasts(Statuspregnant, Statusvirgin) %>%
#'   model_bayes()
model_bayes <- function(.fit, .f = limma::eBayes, ...) {
  .args <- rlang::enexprs(...)
  
  rlang::eval_tidy(rlang::expr(.f(fit = .fit,
                                  !!! .args)))
}

#' Tidying methods for MArrayLM objects
#'
#' These methods tidy the results of differential expression analysis objects.
#' Currently only `limma` is supported.
#'
#' @param x a differential expression fit object. Currently supports `MArrayLM`
#' from the `limma` package.
#' @param conf.int logical. Include confidence intervals?
#' @param exponentiate logical. Should the estimates and (if `conf.int` = `TRUE`)
#' confidence intervals be exponentiated?
#' @param ... additional arguments
#'
#' @return a `tbl`
#'
#' @export
#'
#' @examples
#' counts <- readr::read_delim("data/GSE60450_Lactation-GenewiseCounts.txt", delim = "\t")
#' meta <- readr::read_delim("data/SampleInfo_Corrected.txt", delim = "\t") %>%
#'   mutate(FileName = stringr::str_replace(FileName, "\\.", "-"))
#'
#' id <- as.character(counts$EntrezGeneID)
#' out <- check_sample_names(counts, c(1,2), meta, FileName) %>%
#'   purrr::pluck("mod_count") %>%
#'   filter_genes(., id, "edgeR") %>%
#'   make_voom(., my_design) %>%
#'   model_limma() %>%
#'   make_contrasts(Statuspregnant, Statusvirgin) %>%
#'   model_bayes()
#'
#' out %>% tidy.marray.lm()
#'
#' @rdname tidiers
tidy.marray.lm <- function(x, conf.int = TRUE, exponentiate = FALSE, ...) {
  if (!inherits(x, "MArrayLM")) stop("`x` must be of class `MArrayLM`")
  
  x <- purrr::map(x, ~as.vector(.x))
  margin_error <- sqrt(x$s2.post)*
    x$stdev.unscaled*qt(0.975, df = x$df.total)
  results <- tibble(gene = unlist(x$genes),
                    estimate = x$coefficients,
                    std.error = x$s2.post,
                    statistic = x$t,
                    p.value = x$p.value,
                    conf.low = x$coefficients - margin_error,
                    conf.high = x$coefficients + margin_error,
                    stringsAsFactors = FALSE)
  
  if (exponentiate) {
    results$estimate <- exp(results$estimate)
    results$conf.low <- exp(results$conf.low)
    results$conf.high <- exp(results$conf.high)
  }
  
  if (!conf.int) {
    results <- dplyr::select(results, -conf.low, -conf.high)
  }
  
  results
}
