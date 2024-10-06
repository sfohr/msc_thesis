setwd("~/projects/thesis_dev/BacSC/analysis/P_aero_S2S3")
library(dplyr)
library(ggplot2)
library(latex2exp)
library(forcats)
library(cowplot)

# el <- read.delim("../../data/control_expression_classification.tsv")
path <- "./theta_3BNMD_r10_1000iters_raw_counts_long_expr_lvl_classes.csv"

# visualize annotated anmd reconstruction to see if it can serve for imputation
d <- readr::read_csv(path) |> 
  rename("id"="...1") |> 
  mutate(
#    expression_level = case_when(is.na(expression_level) ~ "not_expressed", .default=expression_level),
    expression_level = factor(expression_level, 
                              levels = c("not_expressed", "low", "medium_low", 
                                         "medium_high", "high"), 
                              labels = c("not expressed", "low", "medium low", 
                                         "medium high", "high"),
                              ordered = TRUE)
  )

table(d$expression_level, useNA = "always")

# missing loci
table(d$locus_tag[is.na(d$expression_level)])


expr_class <- d |> 
  select(locus_tag, expression_level) |> 
  distinct()

result <- result |> 
  left_join(expr_class)

prop.table(table(result$nmd_class, result$expression_level), margin = 1)*100

set.seed(687418221)
#sample_frac <- 0.01
sample_frac <- 0.1
sample_idx <- sample.int(nrow(d), size = as.integer(sample_frac * nrow(d)), replace = FALSE)


d_sub <- d[sample_idx, ]

# d_sub |> 
#   filter(nmd_expression > 0) |> 
#   mutate(nmd_expression = log(nmd_expression)) |> 
#   #mutate(nmd_expression = )
#   ggplot(aes(x=expression_level, y=nmd_expression, color=expression_level)) +
#   geom_jitter(alpha = 0.2) +
#   geom_violin() +
#   scale_fill_brewer(palette = "Dark2")
#   #geom_point(alpha=0.8)

# d_sub <- d_sub |> tidyr::pivot_longer(
#   cols = starts_with("nmd_expression_"), 
#   names_pattern = "nmd_expression_(.*)iter", 
#   names_to = "iterations",
#   values_to = "nmd_expression") |> 
#   mutate(iterations = as.integer(iterations))

#expr_lvl_recode <- c("not_expressed", "low", "medium_low", "medium_high", "high")

expr_lvl_hist <- d |> 
  #filter(iterations == 500) |> 
  select(expression_level, nmd_expression) |> 
  #mutate(expression_level = fct_recode(expression_level, "not expressed" = "not_expressed", "medium low" = "medium_low", "medium high" = "medium_high")) |> 
  tidyr::drop_na() |> 
  ggplot(aes(x=nmd_expression, group=expression_level)) +
  geom_histogram(aes(fill=expression_level), binwidth = 0.1) +
  geom_vline(xintercept = 0, color = "white") +
  coord_cartesian(ylim = c(0, 25)) +
  #coord_cartesian(ylim = c(0, 25)) +
  facet_grid(expression_level~.) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = TeX("$\\Theta_{3B-NMD}$"), y = "Counts") +
  guides(fill = "none") +
  theme(strip.text.y = element_blank())

  # labs(
  #   title = TeX("Distribution of NMD reconstruction $\\Theta$"), 
  #   subtitle = "3B-NMD on raw count expression matrix, Pseudomonas S3 (regular iron env), 500 iters, latent rank 5", 
  #   x = TeX("$\\Theta_{3B-NMD}$"), y = "Counts", fill = "Expression level class")

# zoomed in
expr_lvl_hist_zoomed <- d |> 
  #filter(iterations == 500) |> 
  select(expression_level, nmd_expression) |> 
  # mutate(expression_level = fct_recode(expression_level, "not expressed" = "not_expressed", "medium low" = "medium_low", "medium high" = "medium_high")) |> 
  tidyr::drop_na() |> 
  ggplot(aes(x=nmd_expression, group=expression_level)) +
  geom_histogram(aes(fill=expression_level), binwidth = 0.1) +
  geom_vline(xintercept = 0, color = "white") +
  coord_cartesian(ylim = c(0, 25), xlim = c(-25, 20)) +
  facet_grid(expression_level~.) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = TeX("$\\Theta_{3B-NMD}$"), y = "", fill = "Expression level\nclass") +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks = element_blank())

expr_lvl_plots <- plot_grid(
  expr_lvl_hist + theme(legend.position="none"), 
  expr_lvl_hist_zoomed + theme(legend.position="none"), 
  labels = "AUTO")

legend <- get_legend(
  # create some space to the left of the legend
  expr_lvl_hist_zoomed + theme(legend.box.margin = margin(0, 0, 0, 2))
)


plot_ready <- plot_grid(expr_lvl_plots, legend, rel_widths = c(2, .4))
plot_ready

ggsave2("./nmd_reconstruction_imputation_hist_3BNMD_r10_04_counts.pdf", 
        plot = plot_ready, width = 12, height = 6)
