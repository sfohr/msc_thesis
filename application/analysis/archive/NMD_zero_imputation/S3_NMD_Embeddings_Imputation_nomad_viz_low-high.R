library(ggplot2)
library(dplyr)
library(ggdark)
library(latex2exp)

path <- "./data/gene_pair_df_4_viz_200genes_low-high_withoutRepl.csv"

d <- read.csv(path)

# anzahl expression values kleiner null pro gene pair & gene
neg_val_count_by_gene_pair <- d |> 
  group_by(gene_pair_id, Gene) |> 
  summarise(n_vals_below_zero = sum(Expression < 0.0)) |> 
  ungroup() |> 
  group_by(gene_pair_id) |> 
  summarise(lowest_count = min(n_vals_below_zero))

gene_pairs_bigger_5 <- neg_val_count_by_gene_pair$gene_pair_id[neg_val_count_by_gene_pair$lowest_count > 5]

d_agg <- d |> 
  #
  #filter(gene_pair_id %in% gene_pairs_bigger_5) |> 
  #filter(Expression <= 0) |> 
  group_by(gene_pair_id, Gene) |> 
  summarise(
    median_expression = median(Expression),
    mean_expression = mean(Expression),
    mean_negative_expression = mean(Expression[Expression < 0.0]),
    median_negative_expression = median(Expression[Expression < 0.0]),
    median_negative_expression = tidyr::replace_na(median_negative_expression, 0.0),
    n_negative = sum(Expression < 0.0),
    n = n(),
    share_negative = n_negative / n,
    min_val = min(Expression),
    dist = unique(dist),
    gene_pair_id = unique(gene_pair_id),
    expression_lvl_class = unique(expression_lvl_class)
  )

pair_diffs <- d_agg |> 
  group_by(gene_pair_id) |> 
  arrange(desc(expression_lvl_class)) |> 
  tidyr::pivot_wider(
    id_cols = gene_pair_id, 
    names_from = "expression_lvl_class", 
    values_from = c("median_expression", "mean_expression", "share_negative", "min_val", "n_negative", "mean_negative_expression", "median_negative_expression")) |> 
  mutate(
    median_diff = median_expression_low - median_expression_high,
    mean_diff = mean_expression_low - mean_expression_high,
    
    median_negative_expression_low = tidyr::replace_na(median_negative_expression_low, 0.0),
    median_negative_expression_high = tidyr::replace_na(median_negative_expression_high, 0.0),
    median_neg_diff = median_negative_expression_low - median_negative_expression_high,
    
    mean_negative_expression_low = tidyr::replace_na(mean_negative_expression_low, 0.0),
    mean_negative_expression_high = tidyr::replace_na(mean_negative_expression_high, 0.0),
    mean_neg_diff = mean_negative_expression_low - mean_negative_expression_high,
    share_negative_diff = share_negative_low - share_negative_high,
    n_negative_diff = n_negative_low - n_negative_high) |> 
  ungroup()

pair_diffs |>
  mutate(y_order = rank(mean_diff)) |>
  ggplot(aes(x = mean_diff, y = y_order, color = mean_diff)) + geom_point()

# connected dots: cool
d_agg |>
  ungroup() |>
  left_join(pair_diffs, by="gene_pair_id") |>
  filter(dist == 0) |>
  mutate(y_order = dense_rank(desc(median_neg_diff))) |>
  arrange(desc(dist)) |>
  ggplot(aes(x=mean_negative_expression, y = y_order, group = gene_pair_id,
             color=expression_lvl_class)) +
  geom_line() + geom_point()


# nur die kleinen striche: EXZELLENT!
d_plt <- d |> 
  filter(dist == 0.0) |> 
  left_join(pair_diffs, by="gene_pair_id") |> 
  mutate(
    y_order = dense_rank(desc(mean_neg_diff)),
    y_order = case_when(
      expression_lvl_class == "low" ~ y_order + 0.15,
      expression_lvl_class == "high" ~ y_order - 0.15)
  )

p <- d_plt |> 
  ggplot(aes(x = Expression, y = y_order, color = expression_lvl_class)) +
  geom_point(shape="|", alpha = 0.6) +
  theme(panel.background = element_rect(fill="black"), 
        panel.grid = element_blank(), 
        panel.grid.major.x = element_line(color="darkgrey", linewidth = 0.05),
        panel.grid.minor.x = element_line(color="darkgrey", linewidth = 0.05))
#dark_theme_linedraw()

p +  labs(
  title = TeX("Distribution of NMD reconstruction $\\Theta$ on paired genes"), 
  subtitle = "3B-NMD on raw count expression matrix, Pseudomonas S3 (regular iron env), 500 iters, latent rank 5. \nGene pairs matched by mean and sparsity (exact match).", 
  x = TeX(r'($\Theta_{3B-NMD}$)'), 
  y = TeX(r"(rank$\left( \bar{\Theta^{(-)}_{low}} - \bar{\Theta^{(-)}_{high}}\right)$)"), 
  color = "Expression level class") +
  guides(color = guide_legend( 
    override.aes=list(shape = 15))) +
  geom_vline(xintercept = 0, color="black")

