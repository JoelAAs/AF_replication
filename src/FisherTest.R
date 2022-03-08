library(tidyverse)

#### FUNCtions
row_mean_impute <- function(geno_mat){
  for (i in 2:nrow(geno_mat)) {
    na_idx <- is.na(geno_mat[i, ])
    if (any(na_idx)) {
      imp_value <- mean(as.numeric(geno_mat[i, -1]), na.rm = T) %>% round()
      geno_mat[i, na_idx] <- imp_value
    }
  }
  return(geno_mat)
}

get_driver_variants <- function(geno_mat, response_vector) {
  if (ncol(geno_mat) == 0) {
    return(NA)
  }
  res <- list()
  colnames <- colnames(geno_mat)
  for (i in 1:ncol(geno_mat)) {
    pk <- colnames[i]
    res[[i]] <- format_output(geno_mat[, -i], response_vector)
    res[[i]]["PK"] <- pk
  }
  return(bind_rows(res))
}


fisher_test <- function(geno_mat, response_vector){
  carriers <- as.numeric(apply(geno_mat,1,FUN=function(x) any(x>=1)))
  fisher_res <- fisher.test(carriers, response_vector)
  return(fisher_res)
}

format_output <- function(geno_mat, response_vector){
  if (!is.matrix(geno_mat)){
    geno_mat <- as.matrix(geno_mat)
  }
  
  names <- c(
    "OR.fisher",
    "OR.95CI.lower_fisher",
    "OR.95CI.upper_fisher",
    "pvalue.fisher",
    "n.carrier.case",
    "n.carrier.case.pct",
    "n.carrier.control",
    "n.carrier.control.pct"
    
  )
  
  if (ncol(geno_mat) >= 1) {
    fisher_res <- fisher_test(geno_mat, response_vector)
    n_carrires_cases    <- sum(as.logical(rowSums(geno_mat))[as.logical(response_vector)])
    n_carrires_controls <- sum(as.logical(rowSums(geno_mat))[as.logical(!response_vector)])
    
    results <-  c(
        fisher_res$estimate,
        fisher_res$conf.int[1],
        fisher_res$conf.int[2],
        fisher_res$p.value,
        n_carrires_cases,
        n_carrires_cases/sum(response_vector),
        n_carrires_controls,
        n_carrires_controls/sum(!as.logical(response_vector))
      )
    
  } else {
    results <- rep(NA, length(names))
  }
  names(results) <- names
  return(results)
}




args <- commandArgs(trailingOnly = T)
genotype_matrix_filename  <- args[1]
variant_ann_filename      <- args[2]
phenotype_vector_filename <- args[3]
gene_name                 <- args[4]
output_test               <- args[5]
output_driver             <- args[6]


genotype_matrix      <- read.csv(genotype_matrix_filename, sep = "\t", stringsAsFactors = F, check.names = F)
variant_ann          <- read.csv(variant_ann_filename, sep = "\t", stringsAsFactors = F, check.names = F)
phenotype            <- read.csv(phenotype_vector_filename, sep ="\t", stringsAsFactors = F, check.names = F)

genotype_matrix <- genotype_matrix[, c("PK", phenotype$sample.id)]
genotype_matrix <- row_mean_impute(genotype_matrix)
sample_vf       <- genotype_matrix[, -1] %>% 
  {rowSums(.[, phenotype$sample.id])/length(phenotype$sample.id)}

non_sigular_variants <- genotype_matrix$PK[!sample_vf %in% c(1, 0)]

variant_ann <- variant_ann %>% 
  filter(PK %in% non_sigular_variants)
genotype_matrix <- genotype_matrix %>% 
  filter(PK %in% variant_ann$PK)
geno_mat <- t(genotype_matrix[,-1])

if (!is_empty(geno_mat)){
  colnames(geno_mat) <- genotype_matrix$PK
}

fisher_results <- format_output(geno_mat, phenotype$case) %>%
  as.list() %>% 
  data.frame() %>%
  add_column(gene = gene_name) %>% 
  write.table(output_test, sep="\t", row.names = F)

get_driver_variants(geno_mat, phenotype$case) %>%
  merge(variant_ann, on="PK") %>%
  write.table(output_driver, sep="\t", row.names = F)





