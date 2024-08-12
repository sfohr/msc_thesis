library(dplyr)
library(ggplot2)
library(tidyr)

d <- readRDS("./S3_filtered_with_counts.RDS")

set.seed(64871)
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

df_sub <- d[keep, -1][ , total_counts_per_gene > quantile(total_counts_per_gene, genes_out)]
cl_abs <- readRDS("./cluster_solution_s3_all_absolute_abundance_10perc_k5.RDS")

table(cl_abs$cluster)
df_sub$cluster_abs <- cl_abs$cluster

df_sub %>% 
  mutate(id = seq_len(n())) %>% 
  group_by(cluster_abs) %>% 
  tidyr::pivot_longer(cols = starts_with("PA"), 
                      names_to = "gene", values_to = "counts") %>% 
  group_by(cluster_abs, gene) %>% 
  summarise(avg = mean(counts))

df_top_genes <- df_sub %>% 
  mutate(id = seq_len(n())) %>%
  pivot_longer(cols = starts_with("PA"), names_to = "gene", values_to = "counts") %>%
  group_by(cluster_abs, gene) %>%
  summarise(avg_counts = mean(counts), .groups = 'drop') %>%
  group_by(cluster_abs) %>% 
  mutate(rank = rank(-avg_counts, ties.method = "first")) %>%
  filter(rank <= 10) %>%
  ungroup()

df_top_genes_per_cluster <- df_sub %>%
  pivot_longer(cols = starts_with("PA"), names_to = "gene", 
               values_to = "counts") %>%
  group_by(cluster_abs, gene) %>%
  summarise(avg_counts = mean(counts)) %>%
  #group_by(cluster_abs) %>%
  ungroup(gene) %>% 
  mutate(rank = dense_rank(desc(avg_counts))) %>%
  filter(rank <= 10) %>%
  ungroup() %>%
  arrange(cluster_abs, desc(avg_counts))

ggplot(df_top_genes, 
       aes(x = reorder(gene, -avg_counts), 
           y = avg_counts, fill = factor(cluster_abs))) +
  facet_grid(cluster_abs~.) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(x = "Gene", y = "Average Expression", title = "Top 10 Genes per Cluster by Average Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
