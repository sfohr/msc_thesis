#library(PHclust)
#library(readxl)

setwd("~/ph")
source("./Hybrid.R")
# d <- read_xlsx("./S3_filtered_with_counts.xlsx")
# saveRDS(d, "./S3_filtered_with_counts.RDS")

logfile <- paste0("log_", lubridate::now(), ".txt")

add_date_at_start <- function(msg) {
  now <- lubridate::now()
  return(paste0(now, " -> ", msg))
}

write_log <- function(line, append = TRUE, path = logfile) {
  write(line, file = path, append = append)     
}

write_log(add_date_at_start("Read data"), append = FALSE)

#### Read data ----
d <- readRDS("./S3_filtered_with_counts.RDS")

seed <- 64871
set.seed(seed)
samples_in <- 0.5
samples_out <- 1.0 - samples_in
keep <- sample(c(TRUE, FALSE), size = nrow(d), replace = TRUE, 
               prob = c(samples_in, samples_out))
dm <- as.matrix(d[keep, -1])
total_counts_per_gene <- colSums(dm)
genes_in <- 1.0
genes_out <- 1.0 - genes_in
dm_high <- dm[ , total_counts_per_gene > quantile(total_counts_per_gene, genes_out)]
dim(dm_high)

rows2keep <- which(keep)
cols2keep <- which(total_counts_per_gene > quantile(total_counts_per_gene, genes_out))

k <- 4

#### cluster for absolute abundances ----
write_log(
  add_date_at_start(
    paste0("Do clustering with ", k, " clusters (absolute abundance")), 
  append = TRUE)

clust <- PHcluster(dm_high, Treatment = rep_len(1, ncol(dm_high)), 
                   nK=k, 
                   absolute=TRUE, method = "EM")

#### saving rds for absolute abundances ----
write_log(add_date_at_start("Saving result to file"), 
          append = TRUE)
abs_path <- paste0(
  "cluster_solution_s3_all_absolute_abundance_", 
  samples_in, "samples_", 
  genes_in, "genes_", 
  "k", k,".RDS")
saveRDS(clust, abs_path)

# abs_cl <- readRDS("./cluster_solution_s3_all_absolute_abundance_0.6samples_0.8genes_k5.RDS")
# 
# writeLines(jsonlite::toJSON(list("rows2keep" = rows2keep, 
#                                  "cols2keep" = cols2keep, 
#                                  "cluster_abs" = abs_cl$cluster,
#                                  "seed" = seed)), "./05samples_08genes.json")
# jsonlite::toJSON(list("rows2keep" = rows2keep, "cols2keep" = cols2keep, "seed" = seed))
# 

write_log(
  add_date_at_start(
    paste0("Do clustering with ", k, " clusters (relative abundance")),
  append = TRUE)
#### cluster for relative abundances ----
k_rel <- k
clust_rel <- PHcluster(dm_high, Treatment = rep_len(1, ncol(dm_high)), 
                       nK=k_rel, absolute=FALSE, method = "EM")

#### saving rds for relative abundances ----
rel_path <- paste0(
  "cluster_solution_s3_all_relative_abundance_", 
  samples_in, "samples_", 
  genes_in, "genes_", 
  "k", k_rel,".RDS")
saveRDS(clust_rel, rel_path)