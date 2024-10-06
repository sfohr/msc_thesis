setwd("~/projects/thesis_dev/BacSC/analysis/P_aero_S2S3")
library(ggplot2)
library(dplyr)
library(ggdark)
library(latex2exp)
library(cowplot)


path_rank10 <- "../../data/gene_pair_df_4_viz_100genes_3BNMD_rank10_low-high_mean_and_sparsity_whereXzero.csv"

d_r10 <- read.csv(path_rank10) |> filter(!is.na(Expression))

d_n_neg <- d_r10 |> 
  group_by(gene_pair_id, Gene) |> 
  summarise(
    n_negative = sum(Expression < 0.0)) |> 
  ungroup() |> 
  group_by(gene_pair_id) |> 
  summarise(has_n_smaller_3 = any(n_negative < 3))

ttests <- d_r10 |> 
  mutate(expression_lvl_class = ordered(expression_lvl_class, c("low", "high"), c("low", "high"))) |> 
  left_join(d_n_neg, by="gene_pair_id") |> 
  filter(Expression <= 0, !has_n_smaller_3) |> 
  group_by(gene_pair_id) %>% 
  summarise(t_test = t.test(Expression ~ expression_lvl_class, 
                            alternative = "two.sided") |> broom::tidy()) |> 
  tidyr::unnest(cols = c(t_test))
  # summarise(t_test = t.test(x=Expression[expression_lvl_class == "low"],
  #                           y=Expression[expression_lvl_class == "high"], 
  #                           alternative = "greater") |> broom::tidy())

ttests |> 
  left_join(pair_diffs) |> 
  mutate(y_order = dense_rank(desc(mean_neg_diff)),
         col = estimate < 0) |>
  ggplot(aes(x=estimate, y=y_order)) +
  geom_point(aes(color=col)) +
  geom_errorbar(aes(xmin=conf.low, xmax= conf.high, color=col))


d_agg_r10 <- d_r10 |> 
  group_by(gene_pair_id, Gene) |> 
  summarise(
    median_expression = median(Expression),
    mean_expression = mean(Expression),
    mean_negative_expression = mean(Expression[Expression < 0.0]),
    lci = get_confint(Expression[Expression < 0.0], type="lower"),
    uci = get_confint(Expression[Expression < 0.0], type="upper"),
    mean_negative_expression = tidyr::replace_na(mean_negative_expression, 0.0),
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

pair_diffs <- d_agg_r10 |> 
  group_by(gene_pair_id) |> 
  arrange(desc(expression_lvl_class)) |> 
  tidyr::pivot_wider(
    id_cols = gene_pair_id, 
    names_from = "expression_lvl_class", 
    values_from = c("median_expression", "mean_expression", "lci", "uci", "share_negative", 
                    "min_val", "n_negative", "mean_negative_expression", 
                    "median_negative_expression")) |> 
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

# connected dots
x_lim_negative <- min(d_agg_r10$mean_negative_expression) - 0.02
connect_dotplot_r10_mean <- d_agg_r10 |>
  ungroup() |>
  left_join(pair_diffs, by="gene_pair_id") |>
  #filter(dist == 0) |>
  mutate(y_order = dense_rank(desc(mean_neg_diff))) |>
  arrange(desc(dist)) |>
  ggplot(aes(x=mean_negative_expression, y = y_order, group = gene_pair_id,
             color=expression_lvl_class)) +
  geom_line() + 
  geom_point() + 
  ylab("") +
  xlab(TeX(r"(mean$( \Theta^{(-)})$)")) +
  labs(color = "Expression\nLevel Class") +
  xlim(c(x_lim_negative, 0))


connect_dotplot_r10_mean_cis <- d_agg_r10 |>
  ungroup() |>
  left_join(pair_diffs, by="gene_pair_id") |>
  #filter(dist == 0) |>
  mutate(y_order = dense_rank(desc(mean_neg_diff))) |>
  arrange(desc(dist)) |>
  ggplot(aes(x=mean_negative_expression, y = y_order, group = gene_pair_id,
             color=expression_lvl_class)) +
  #geom_line() + 
  geom_point(aes(x=mean_negative_expression, y=y_order, color=expression_lvl_class)) +
  geom_errorbar(aes(xmin=lci, xmax= uci, group=Gene, color=expression_lvl_class)) +
  ylab("") +
  xlab(TeX(r"(mean$( \Theta^{(-)})$)")) +
  labs(color = "Expression\nLevel Class")
  #xlim(c(x_lim_negative, 0))



#ggsave2("./visualize_nmd_imputation/low_high_3BNMD_r10_04_connected_dots_mean_diffs.svg",
#        connect_dotplot_r10_mean)

connect_dotplot_r10_mean

connect_dotplot_r10_median <- d_agg_r10 |>
  ungroup() |>
  left_join(pair_diffs, by="gene_pair_id") |>
  #filter(dist < 0.04) |>
  mutate(y_order = dense_rank(desc(mean_neg_diff))) |>
  arrange(desc(dist)) |>
  ggplot(aes(x=median_negative_expression, y = y_order, group = gene_pair_id,
             color=expression_lvl_class)) +
  geom_line() + 
  geom_point() + 
  ylab(TeX(r"(rank$\left( \bar{\Theta}^{(-)}_{low} - \bar{\Theta}^{(-)}_{high}\right)$)")) +
  xlab(TeX(r"(median$( \Theta^{(-)})$)")) +
  xlim(c(x_lim_negative, 0))


connect_dotplot_r10_median
#ggsave2("./visualize_nmd_imputation/low_high_3BNMD_r10_04_connected_dots_mean_diffs.svg")

# combine median and mean plot
connected_dot_plots <- plot_grid(
  connect_dotplot_r10_median + theme(legend.position="none"), 
  connect_dotplot_r10_mean + theme(legend.position="none"), 
labels = "AUTO")

legend <- get_legend(
  # create some space to the left of the legend
  connect_dotplot_r10_mean + theme(legend.box.margin = margin(0, 0, 0, 1))
)


connected_dots_plot_ready <- plot_grid(connected_dot_plots, legend, rel_widths = c(2, 0.3))
connected_dots_plot_ready

ggsave2("./visualize_nmd_imputation/low_high_3BNMD_r10_04_connected_dots_multiplot_diffs.pdf", 
        plot = connected_dots_plot_ready, width = 10, height = 5)


# fance kleine striche plot ----
d_plt_lines <- d_r10 |> 
  #filter(dist == 0.0, Expression < 0.0) |> 
  left_join(pair_diffs, by="gene_pair_id") |> 
  mutate(
    y_order = dense_rank(desc(mean_neg_diff)),
    y_order = case_when(
      # expression_lvl_class == "low" ~ y_order + 0.15,
      # expression_lvl_class == "high" ~ y_order - 0.15)
      expression_lvl_class == "low" ~ y_order + 0.25,
      expression_lvl_class == "high" ~ y_order - 0.25)
    )

p <- d_plt_lines |> 
  ggplot(aes(x = Expression, y = y_order, color = expression_lvl_class)) +
  geom_point(shape="|", alpha = 0.6) +
  theme(panel.background = element_rect(fill="black"), 
        panel.grid = element_blank(), 
        panel.grid.major.x = element_line(color="darkgrey", linewidth = 0.05),
        panel.grid.minor.x = element_line(color="darkgrey", linewidth = 0.05))
  #dark_theme_linedraw()

striche <- p +  labs(
  #title = TeX("Distribution of NMD reconstruction $\\Theta$ on paired genes"), 
  #subtitle = "3B-NMD on raw count expression matrix, Pseudomonas S3 (regular iron env), 500 iters, latent rank 5. \nGene pairs matched by mean and sparsity (exact match).", 
  x = TeX(r'($\Theta_{3B-NMD}$)'), 
  y = TeX(r"(rank$\left( \bar{\Theta}^{(-)}_{low} - \bar{\Theta}^{(-)}_{high}\right)$)"), 
  color = "Expression level class") +
  guides(color = guide_legend( 
    override.aes=list(shape = 15))) +
  geom_vline(xintercept = 0, color="black") +
  xlim(c(-15, 5))

cowplot::ggsave2(plot=striche, filename = "./fancy_plot.png", dpi = 300, width = 6, height = 8)
