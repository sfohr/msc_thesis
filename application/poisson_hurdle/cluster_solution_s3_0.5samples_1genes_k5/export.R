# save rows2keep and clustering solution for python

setwd("~/ph/cluster_solution_s3_0.5samples_1genes_k5/")

data_path <- "../S3_filtered_with_counts.RDS"
cl_abs_path <- "./cluster_solution_s3_all_absolute_abundance_0.5samples_1genes_k5.RDS"
cl_rel_path <- "./cluster_solution_s3_all_relative_abundance_0.5samples_1genes_k5.RDS"

d <- readRDS(data_path)
cl_abs <- readRDS(cl_abs_path)
cl_rel <- readRDS(cl_rel_path)

seed <- 64871

make_export <- function(data_path, 
                        cluster_path, 
                        seed, 
                        sample_frac=1.0, 
                        gene_frac=1.0, 
                        absolute = TRUE) {
  d <- readRDS(data_path)
  res <- readRDS(cluster_path)
  
  set.seed(seed)
  samples_out <- 1.0 - sample_frac
  keep <- sample(c(TRUE, FALSE), size = nrow(d), replace = TRUE, 
                 prob = c(sample_frac, samples_out))
  dm <- as.matrix(d[keep, -1])
  total_counts_per_gene <- colSums(dm)
  genes_out <- 1.0 - gene_frac
  dm_high <- dm[ , total_counts_per_gene > quantile(total_counts_per_gene, genes_out)]
  
  rows2keep <- which(keep)
  cols2keep <- which(total_counts_per_gene > quantile(total_counts_per_gene, genes_out))
  
  res$rows2keep <- rows2keep
  res$cols2keep <- cols2keep
  res$absolute <- absolute
  return(res)
}

exp_abs <- make_export(data_path, 
                   cluster_path = cl_abs_path, seed = seed, 
                   sample_frac = 0.5, gene_frac = 1.0, absolute = TRUE)

exp_rel <- make_export(data_path, 
                       cluster_path = cl_rel_path, seed = seed, 
                       sample_frac = 0.5, gene_frac = 1.0, absolute = FALSE)

jsonlite::write_json(exp_abs, "./S3_abs_05samples_1genes_k5.json")
jsonlite::write_json(exp_rel, "./S3_rel_05samples_1genes_k5.json")
