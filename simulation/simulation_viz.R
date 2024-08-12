library(ggplot2)
library(dplyr)
library(latex2exp)
library(cowplot)

d <- readRDS("./all_sims_lrz.RDS")

d <- d |> 
  filter(m != 5000) |> 
  mutate(algorithm_factor = factor(algorithm, c("A_NMD", "3B_NMD"), labels = c("A-NMD", "3B-NMD")))

filter_by_epsilon_convergence <- function(df, niters = 10, epsilon = 1e-4) {
  d |>
    group_by(m, r, sparsity_level, replication, initialization) |>
    mutate(
      has_converged_epsilon = abs(
        rel_errors - lag(rel_errors, niters)
      ) < epsilon
    ) |>
    ungroup() |>
    mutate(
      has_converged_epsilon = tidyr::replace_na(has_converged_epsilon, FALSE)
    ) |>
    filter(!has_converged_epsilon)
}

rel_error_threshold <- 0.03

# filter all
d <- d |>
  group_by(identifier) |> 
  mutate(converged_by_error = cumsum(cumsum(rel_errors < rel_error_threshold))) |> 
  filter(converged_by_error < 2) |> 
  ungroup()

# summary data frame
d_sum <- d |> 
  group_by(identifier) |> 
  arrange(desc(-iterations)) |>
  summarise(
    algorithm = unique(algorithm),
    algorithm_factor = unique(algorithm_factor),
    m = unique(m),
    r = unique(r),
    sparsity_level = unique(sparsity_level),
    actual_sparsity = unique(actual_sparsity),
    replication = unique(replication),
    initialization = unique(initialization),
    final_rel_error = last(rel_errors),
    total_time = sum(times),
    mean_time_per_iter = mean(times),
    median_time_per_iter = median(times),
    n_iters = max(iterations),
    converged_by_rel_error = max(converged_by_error)
  ) |> 
  ungroup()

#### DIMENSION ----

# total log10 time until convergence over dimension
dimension_X_totaltime <- d_sum |> 
  ggplot(aes(x = factor(m), y = total_time, fill = algorithm_factor)) +
  geom_violin(draw_quantiles = c(0.5)) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Total time (in secs)", x = "Matrix dimension", fill = "Algorithm") +


ggsave("dimension_X_totaltime.pdf", 
       plot=dimension_X_totaltime, device = "pdf", width = 7, height = 5, 
       units = "in", dpi = 300, path = "./plots/")

summary(lm(log(total_time) ~ m + r + actual_sparsity*algorithm, data=d_sum))
summary(lm(log(mean_time_per_iter) ~ m + r + algorithm*actual_sparsity, data=d_sum))
qplot(predict(lm(log(total_time) ~ m + algorithm, data=d_sum)), log(d_sum$total_time))
# plot(lm(log(total_time) ~ m*r*actual_sparsity, data=d_sum))

# total log10 iterations until convergence over dimension
dimension_X_niterations <- d_sum |> 
  ggplot(aes(x = factor(m), y = n_iters, fill = algorithm_factor)) +
  geom_violin(draw_quantiles = c(0.5)) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Iterations until convergence", x = "Matrix dimension", fill = "Algorithm")

ggsave("dimension_X_niterations.pdf", 
       plot=dimension_X_niterations, device = "pdf", width = 7, height = 5, 
       units = "in", dpi = 300, path = "./plots/")

# total log10 time per iteration over dimension
# timeperiter_X_dimension <- d_sum |> 
#   ggplot(aes(x = factor(m), y = mean_time_per_iter, fill = algorithm_factor)) +
#   geom_boxplot() +
#   #scale_y_log10(labels = function(secs) round(secs, 2)) +
#   labs(y = "Time per iteration (in seconds)", x = "Matrix dimension", fill = "Algorithm")

dimension_X_timeperiter <- d |> 
  filter(times > 0) |> 
  ggplot(aes(x = factor(m), y = times, fill = algorithm_factor)) +
  geom_violin(draw_quantiles = c(0.5)) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Time per iteration (in secs)", x = "Matrix dimension", fill = "Algorithm")

ggsave("dimension_X_timeperiter.pdf", 
       plot=dimension_X_timeperiter, device = "pdf", width = 7, height = 5, 
       units = "in", dpi = 300, path = "./plots/")

dimension_multiplot <- ggpubr::ggarrange(
  dimension_X_niterations, dimension_X_totaltime, dimension_X_timeperiter,
  align = "h", ncol = 3, labels = c("A", "B", "C"),
  common.legend = TRUE, legend = "bottom"
)

ggsave("dimension_multiplot.pdf", 
       plot=dimension_multiplot, device = "pdf", width = 10, height = 3.5, 
       units = "in", dpi = 300, path = "./plots/")


#### LATENT RANK ----
rank_X_totaltime <- d_sum |> 
  ggplot(aes(x = factor(r), y = total_time, fill = algorithm_factor)) +
  geom_violin(draw_quantiles = c(0.5)) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Total time (in secs)", x = "Rank", fill = "Algorithm")

ggsave("rank_X_totaltime.pdf", 
       plot=rank_X_totaltime, device = "pdf", width = 7, height = 5, 
       units = "in", dpi = 300, path = "./plots/")

summary(lm(log(total_time) ~ m*r + actual_sparsity, data=d_sum))
# plot(lm(log(total_time) ~ m*r*actual_sparsity, data=d_sum))

# total log10 iterations until convergence over dimension
rank_X_niterations <- d_sum |> 
  ggplot(aes(x = factor(r), y = n_iters, fill = algorithm_factor)) +
  geom_violin(draw_quantiles = c(0.5)) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Iterations until convergence", x = "Rank", fill = "Algorithm")

ggsave("rank_X_niterations.pdf", 
       plot=rank_X_niterations, device = "pdf", width = 7, height = 5, 
       units = "in", dpi = 300, path = "./plots/")

rank_X_timeperiter <- d |> 
  filter(times > 0) |> 
  ggplot(aes(x = factor(r), y = times, fill = algorithm_factor)) +
  geom_violin(draw_quantiles = c(0.5)) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Time per iteration (in secs)", x = "Rank", fill = "Algorithm")

ggsave("rank_X_timeperiter.pdf", 
       plot=rank_X_timeperiter, device = "pdf", width = 7, height = 5, 
       units = "in", dpi = 300, path = "./plots/")

rank_multiplot <- ggpubr::ggarrange(
  rank_X_niterations, rank_X_totaltime, rank_X_timeperiter,
  align = "h", ncol = 3, labels = c("A", "B", "C"),
  common.legend = TRUE, legend = "bottom"
)

ggsave("rank_multiplot.pdf", 
       plot=rank_multiplot, device = "pdf", width = 10, height = 3.5,  
       units = "in", dpi = 300, path = "./plots/")

#### SPARSITY ----
sparsity_X_totaltime <- d_sum |> 
  ggplot(aes(x = factor(sparsity_level), y = total_time, fill = algorithm_factor)) +
  geom_violin(draw_quantiles = c(0.5)) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Total time (in secs)", x = "Sparsity", fill = "Algorithm")

ggsave("sparsity_X_totaltime.pdf", 
       plot=sparsity_X_totaltime, device = "pdf", width = 7, height = 5, 
       units = "in", dpi = 300, path = "./plots/")

summary(lm(log(total_time) ~ m*r + actual_sparsity, data=d_sum))
# plot(lm(log(total_time) ~ m*r*actual_sparsity, data=d_sum))

# total log10 iterations until convergence over dimension
sparsity_X_niterations <- d_sum |> 
  ggplot(aes(x = factor(sparsity_level), y = n_iters, fill = algorithm_factor)) +
  #geom_violin(draw_quantiles = c(0.5)) +
  geom_boxplot() +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Iterations until convergence", x = "Sparsity", fill = "Algorithm")

ggsave("sparsity_X_niterations.pdf", 
       plot=sparsity_X_niterations, device = "pdf", width = 7, height = 5, 
       units = "in", dpi = 300, path = "./plots/")

sparsity_X_timeperiter <- d |> 
  filter(times > 0) |> 
  ggplot(aes(x = factor(sparsity_level), y = times, fill = algorithm_factor)) +
  geom_boxplot() +
  #geom_violin(draw_quantiles = c(0.5)) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Time per iteration (in secs)", x = "Sparsity", fill = "Algorithm")

ggsave("sparsity_X_timeperiter.pdf", 
       plot=sparsity_X_timeperiter, device = "pdf", width = 7, height = 5, 
       units = "in", dpi = 300, path = "./plots/")

sparsity_multiplot <- ggpubr::ggarrange(
  sparsity_X_niterations, sparsity_X_totaltime, sparsity_X_timeperiter,
  align = "h", ncol = 3, labels = c("A", "B", "C"),
  common.legend = TRUE, legend = "bottom"
)

ggsave("sparsity_multiplot.pdf", 
       plot=sparsity_multiplot, device = "pdf", width = 10, height = 3.5, 
       units = "in", dpi = 300, path = "./plots/")


### algo runs in high dimensions
df_check_oscillations <- d |> 
  group_by(identifier) |> 
  #arrange(iterations) |> 
  mutate(error_increasing = rel_errors < lead(rel_errors)) |> ungroup()
  ggplot(aes(x = iterations, y = rel_errors, group = identifier, color = factor(sparsity_level))) +
  geom_line() +
  facet_grid(m ~ r)


#### MULTIVARIATE VIEW ----

# log total time by rank x dimension x sparsity x algorithm
dimension_X_sparsity_X_rank_totaltime <- d_sum |> 
  ggplot(aes(x = factor(m), y = total_time, color = algorithm_factor)) + 
  geom_jitter(height = 0, alpha = 0.8) + 
  facet_grid(sparsity_level~r) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Total time (in secs)", x = "Matrix dimension", color = "Algorithm") +
  theme(legend.position="bottom")

plot2save <- gridExtra::grid.arrange(
  dimension_X_sparsity_X_rank_totaltime, 
  top = "Rank", right = "Sparsity")

ggsave("dimension_X_sparsity_X_rank_totaltime.pdf", 
       plot=plot2save, device = "pdf", width = 7, height = 7, 
       units = "in", dpi = 300, path = "./plots/")

# n iterations
dimension_X_sparsity_X_rank_niters <- d_sum |> 
  ggplot(aes(x = factor(m), y = n_iters, color = algorithm_factor)) + 
  geom_jitter(height = 0, alpha = 0.8) + 
  facet_grid(sparsity_level~r) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Iterations until convergence", x = "Matrix dimension", color = "Algorithm") +
  theme(legend.position="bottom")

plot2save <- gridExtra::grid.arrange(dimension_X_sparsity_X_rank_niters, top = "Rank", right = "Sparsity")
ggsave("dimension_X_sparsity_X_rank_niters.pdf", 
       plot=plot2save, device = "pdf", width = 7, height = 7, 
       units = "in", dpi = 300, path = "./plots/")

# time per iter
dimension_X_sparsity_X_rank_timeperiter <- d |> 
  filter(times > 0) |> 
  ggplot(aes(x = factor(m), y = times, fill = algorithm_factor, color = algorithm_factor)) + 
  #geom_jitter(height = 0, alpha = 0.05) + 
  geom_violin(draw_quantiles = c(0.5)) +
  facet_grid(sparsity_level~r) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Time per iteration (in secs)", 
       x = "Matrix dimension", fill = "Algorithm") +
  theme(legend.position="bottom") +
  guides(color = "none")

dimension_X_sparsity_X_rank_timeperiter_grid_annot <- gridExtra::grid.arrange(
  dimension_X_sparsity_X_rank_timeperiter, 
  top = "Rank", right = "Sparsity")

ggsave("dimension_X_sparsity_X_rank_timeperiter.pdf", 
       plot=dimension_X_sparsity_X_rank_timeperiter_grid_annot, 
       device = "pdf", width = 7, height = 7, 
       units = "in", dpi = 300, path = "./plots/")

dimension_X_sparsity_X_rank_timeperiter_sparsity_on_x <- d |> 
  filter(times > 0) |> 
  ggplot(aes(x = factor(sparsity_level), y = times, 
             fill = algorithm_factor, color = algorithm_factor)) + 
  #geom_jitter(height = 0, alpha = 0.05) + 
  geom_violin(draw_quantiles = c(0.5)) +
  facet_grid(m~r) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Time per iteration (in secs)", 
       x = "Sparsity", fill = "Algorithm") +
  theme(legend.position="bottom") +
  guides(color = "none")

dimension_X_sparsity_X_rank_timeperiter_sparsity_on_x_grid_annot <- gridExtra::grid.arrange(
  dimension_X_sparsity_X_rank_timeperiter_sparsity_on_x, 
  top = "Rank", right = "Matrix dimension")

ggsave("dimension_X_sparsity_X_rank_timeperiter_sparsity_on_x.pdf", 
       plot=dimension_X_sparsity_X_rank_timeperiter_sparsity_on_x_grid_annot, 
       device = "pdf", width = 7, height = 7, 
       units = "in", dpi = 300, path = "./plots/")

# multiplot total time and niterations: dimension on x ----
dimension_X_sparsity_X_rank_totaltime <- d_sum |> 
  ggplot(aes(x = factor(m), y = total_time, color = algorithm_factor)) + 
  geom_jitter(height = 0, alpha = 0.8) + 
  facet_grid(sparsity_level~r) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Total time (in secs)", x = "Matrix dimension", color = "Algorithm") +
  theme(legend.position="bottom")

dimension_X_sparsity_X_rank_totaltime_with_grid_annot <- gridExtra::grid.arrange(
  dimension_X_sparsity_X_rank_totaltime, 
  top = "Rank")

dimension_X_sparsity_X_rank_niters <- d_sum |> 
  ggplot(aes(x = factor(m), y = n_iters, color = algorithm_factor)) + 
  geom_jitter(height = 0, alpha = 0.8) + 
  facet_grid(sparsity_level~r) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Iterations until convergence", x = "Matrix dimension", color = "Algorithm") +
  theme(legend.position="bottom")

dimension_X_sparsity_X_rank_niters_with_grid_annot <- gridExtra::grid.arrange(
  dimension_X_sparsity_X_rank_niters, 
  top = "Rank", right = "Sparsity")

dimension_X_sparsity_X_rank_time_niters_multiplot <- 
  plot_grid(dimension_X_sparsity_X_rank_totaltime_with_grid_annot, 
            dimension_X_sparsity_X_rank_niters_with_grid_annot, 
            ncol = 2, labels = c("A", "B"))

ggsave("dimension_X_sparsity_X_rank_time_niters_multiplot.pdf", 
       plot=dimension_X_sparsity_X_rank_time_niters_multiplot, 
       device = "pdf", width = 12, height = 6, 
       units = "in", dpi = 300, path = "./plots/")


## different views for logtime and iter multiplot: rank on x ----
dimension_X_sparsity_X_rank_logtotaltime_xrank <- d_sum |> 
  ggplot(aes(x = factor(r), y = total_time, color = algorithm_factor)) + 
  geom_jitter(height = 0, alpha = 0.8) + 
  facet_grid(sparsity_level~m) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Total time (in secs)", x = "Rank", color = "Algorithm") +
  theme(legend.position="bottom")

dimension_X_sparsity_X_rank_logtotaltime_xrank_with_grid_annot <- gridExtra::grid.arrange(
  dimension_X_sparsity_X_rank_logtotaltime_xrank, 
  top = "Matrix dimension")

dimension_X_sparsity_X_rank_logniters_xrank <- d_sum |> 
  ggplot(aes(x = factor(r), y = n_iters, color = algorithm_factor)) + 
  geom_jitter(height = 0, alpha = 0.8) + 
  facet_grid(sparsity_level~m) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(y = "Iterations until convergence", x = "Rank", color = "Algorithm") +
  theme(legend.position="bottom")

dimension_X_sparsity_X_rank_logniters_xrank_with_grid_annot <- gridExtra::grid.arrange(
  dimension_X_sparsity_X_rank_logniters_xrank, 
  top = "Matrix dimension", right = "Sparsity")

dimension_X_sparsity_X_rank_time_niters_multiplot <- 
  plot_grid(dimension_X_sparsity_X_rank_logtotaltime_xrank_with_grid_annot, 
            dimension_X_sparsity_X_rank_logniters_xrank_with_grid_annot, 
            ncol = 2, labels = c("A", "B"))

ggsave("dimension_X_sparsity_X_rank_time_niters_multiplot_rank_onx.pdf", 
       plot=dimension_X_sparsity_X_rank_time_niters_multiplot, 
       device = "pdf", width = 12, height = 6, 
       units = "in", dpi = 300, path = "./plots/")

#### ANOMALY ----
#### 

# ANMD anomaly ----
anmd_total_time_anomaly <- d_sum |>
  filter(algorithm == "A_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = sparsity_level, y = total_time,
    color = sparsity_level, group = identifier
  )) + 
  geom_jitter(alpha = .5, height = 0, show.legend = FALSE) +
  facet_grid(rows = vars(r), cols = vars(m))  +
  labs(
    x = "Sparsity",
    y = "Total time (in secs)",
    color = "Sparsity"
  )

anmd_total_time_anomaly_with_grid_annot <- gridExtra::grid.arrange(
  anmd_total_time_anomaly, top = "Dimension")


anmd_total_logtime_anomaly <- d_sum |>
  filter(algorithm == "A_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = sparsity_level, y = total_time,
    color = sparsity_level, group = identifier
  )) + 
  geom_jitter(alpha = .5, height = 0, show.legend = FALSE) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  labs(
    x = "Sparsity",
    y = element_blank(),
    color = "Sparsity"
  )

anmd_total_logtime_anomaly_with_grid_annot <- gridExtra::grid.arrange(
  anmd_total_logtime_anomaly, top = "Dimension", right = "Rank")

anmd_anomaly_multiplot <- 
  plot_grid(anmd_total_time_anomaly_with_grid_annot, 
            anmd_total_logtime_anomaly_with_grid_annot, 
            ncol = 2, labels = c("A", "B"))

ggsave("anmd_anomaly_multiplot.pdf", 
       plot=anmd_anomaly_multiplot, device = "pdf", width = 12, height = 6, 
       units = "in", dpi = 300, path = "./plots/")

# 3B-NMD anomaly ----
bnmd_time_till_conv <- d_sum |>
  filter(algorithm == "3B_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = sparsity_level, y = total_time,
    color = sparsity_level, group = identifier
  )) + geom_jitter(alpha = .5, height = 0, show.legend = FALSE) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  labs(
    x = "Sparsity",
    y = "Total time (in secs)",
    color = "Sparsity"
  )

bnmd_total_time_anomaly_with_grid_annot <- gridExtra::grid.arrange(
  bnmd_time_till_conv, top = "Dimension")

bnmd_logtime_till_conv <- d_sum |>
  filter(algorithm == "3B_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = sparsity_level, y = total_time,
    color = sparsity_level, group = identifier
  )) + geom_jitter(alpha = .5, height = 0, show.legend = FALSE) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(
    x = "Sparsity",
    y = element_blank(),
    color = "Sparsity"
  )

bnmd_total_logtime_anomaly_with_grid_annot <- gridExtra::grid.arrange(
  bnmd_logtime_till_conv, top = "Dimension", right = "Rank")


bnmd_anomaly_multiplot <- 
  plot_grid(bnmd_total_time_anomaly_with_grid_annot, 
            bnmd_total_logtime_anomaly_with_grid_annot, 
            ncol = 2, labels = c("A", "B"))

ggsave("3bnmd_anomaly_multiplot.pdf", 
       plot=bnmd_anomaly_multiplot, device = "pdf", width = 12, height = 6, 
       units = "in", dpi = 300, path = "./plots/")

# anmd and 3b logtime anomaly mutiplot ----
anmd_total_logtime_anomaly_anmd_3b <- d_sum |>
  filter(algorithm == "A_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = sparsity_level, y = total_time,
    color = sparsity_level, group = identifier
  )) + 
  geom_jitter(alpha = .5, height = 0, show.legend = FALSE) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  labs(
    x = "Sparsity",
    y = "Total time (in secs)",
    color = "Sparsity"
  )

anmd_total_logtime_anomaly_anmd_3b_with_grid_annot <- gridExtra::grid.arrange(
  anmd_total_logtime_anomaly_anmd_3b, top = "Dimension")


anmd_bnmd_logtime_anomaly_multiplot <- 
  plot_grid(anmd_total_logtime_anomaly_anmd_3b_with_grid_annot, 
            bnmd_total_logtime_anomaly_with_grid_annot, 
            ncol = 2, labels = c("A", "B"))

ggsave("anmd_bnmd_logtime_anomaly_multiplot.pdf", 
       plot=anmd_bnmd_logtime_anomaly_multiplot, device = "pdf", width = 12, height = 6, 
       units = "in", dpi = 300, path = "./plots/")

# anmd and 3b time anomaly mutiplot ----
anmd_total_time_anomaly_anmd_3b <- d_sum |>
  filter(algorithm == "A_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = sparsity_level, y = total_time,
    color = sparsity_level, group = identifier
  )) + 
  geom_jitter(alpha = .5, height = 0, show.legend = FALSE) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  labs(
    x = "Sparsity",
    y = "Total time (in secs)",
    color = "Sparsity"
  )

anmd_total_time_anomaly_anmd_3b_with_grid_annot <- gridExtra::grid.arrange(
  anmd_total_time_anomaly_anmd_3b, top = "Dimension")

bnmd_total_time_anomaly_anmd_3b <- d_sum |>
  filter(algorithm == "3B_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = sparsity_level, y = total_time,
    color = sparsity_level, group = identifier
  )) + geom_jitter(alpha = .5, height = 0, show.legend = FALSE) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  labs(
    x = "Sparsity",
    y = element_blank(),
    color = "Sparsity"
  )

bnmd_total_time_anomaly_with_grid_annot <- gridExtra::grid.arrange(
  bnmd_total_time_anomaly_anmd_3b, top = "Dimension", right = "Rank")

anmd_bnmd_time_anomaly_multiplot <- 
  plot_grid(anmd_total_time_anomaly_anmd_3b_with_grid_annot, 
            bnmd_total_time_anomaly_with_grid_annot, 
            ncol = 2, labels = c("A", "B"))

ggsave("anmd_bnmd_time_anomaly_multiplot.pdf", 
       plot=anmd_bnmd_time_anomaly_multiplot, device = "pdf", width = 12, height = 6, 
       units = "in", dpi = 300, path = "./plots/")




#### TODO ----
# Count the runs per
d |>
  select(m, r, sparsity_level, replication, initialization) |>
  distinct() |>
  select(m, r, sparsity_level) |>
  group_by(m, r, sparsity_level) |>
  summarise(n_runs = n()) |> View()

# boxplot: sparsity vs iterations fill algorithm
d_sum |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(x = sparsity_level, y = n_iters, fill = algorithm)) +
  geom_boxplot() +
  labs(
    subtitle = paste0(
      "Convergence criteria:",
      "relative error < 1e-4, epsilon over 10 iterations < 1e-5,",
      "< 10 000 iterations"
    ),
    x = "Sparsity level",
    y = "Iterations to convergence",
    fill = "Algorithm"
  )

d |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  group_by(identifier, algorithm, sparsity_level) |>
  summarise(total_time = sum(times)) |>
  ggplot(aes(x = sparsity_level, y = total_time, fill = algorithm)) +
  geom_boxplot(varwidth = TRUE) +
  labs(
    subtitle = paste0(
      "Convergence criteria:",
      "relative error < 1e-4, epsilon over 10 iterations < 1e-5,",
      "< 10 000 iterations"
    ),
    x = "Sparsity level",
    y = "Seconds till convergence",
    fill = "Algorithm"
  )

d |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  group_by(identifier, algorithm, sparsity_level) |>
  #summarise(total_time = sum(times)) |>
  summarise(iterations = max(iterations)) |> 
  ggplot(aes(x = sparsity_level, y = iterations, fill = algorithm)) +
  geom_boxplot(varwidth = TRUE) +
  labs(
    subtitle = paste0(
      "Convergence criteria:",
      "relative error < 1e-4, epsilon over 10 iterations < 1e-5,",
      "< 10 000 iterations"
    ),
    x = "Sparsity level",
    y = "Iterations till convergence",
    fill = "Algorithm"
  )



d |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  group_by(identifier, algorithm, sparsity_level, m, r) |>
  summarise(total_time = sum(times)) |>
  ggplot(aes(x = sparsity_level, y = log(total_time), color=algorithm)) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  geom_jitter(alpha=.5, height = 0) +
  geom_smooth() +
  labs(subtitle = paste0(
    "Convergence criteria:",
    "relative error < 1e-4, epsilon over 10 iterations < 1e-5,",
    "< 10 000 iterations"
  ),
    x = "Sparsity level", y = "log(total time in secs)", fill="Algorithm")



# n iterations till convergence
d_sum |>
  filter(algorithm == "A_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = sparsity_level, y = n_iters,
    color = sparsity_level, group = identifier
  )) + geom_jitter(alpha = .5, height = 0, show.legend = FALSE) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  labs(
    x = "Sparsity",
    y = "Iterations",
    color = "Sparsity"
  )

d_sum |>
  filter(algorithm == "3B_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = sparsity_level, y = n_iters,
    color = sparsity_level, group = identifier
  )) + geom_jitter(alpha = .5, height = 0, show.legend = FALSE) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  labs(
    x = "Sparsity",
    y = "Iterations",
    color = "Sparsity"
  ) +
  scale_y_log10()

# total time till convergence


anmd_time_till_conv_plots <- plot_grid(
  anmd_time_till_conv, 
  anmd_logtime_till_conv, 
  labels = "AUTO")

ggsave("anmd_time_till_conv.pdf", 
       plot=anmd_time_till_conv_plots, device = "pdf", width = 10, height = 5, 
       units = "in", dpi = 300, path = "./plots/")



# ANMD whole range
anmd_rel_error_total <- d |>
  filter(algorithm == "A_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = iterations, y = rel_errors,
    color = sparsity_level, group = identifier
  )) +
  geom_line(alpha = .2) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  # xlim(c(0,  50)) +
  # ylim(c(0,  1.25)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(
    # subtitle = "Unequal distribution of runs & higher sparsity runs completely missing for m>=2000",
    x = "Iterations",
    y = TeX("Relative error: $\\|max\\left(0, \\Theta \\right) - X\\|_{F} / \\|X\\|_{F}$"),
    color = "Sparsity"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . , name = "Latent rank", breaks = NULL, labels = NULL)) +
  scale_x_continuous(limits = c(0, 3000), sec.axis = sec_axis(~ . , name = "Dimension", breaks = NULL, labels = NULL))

# 3B whole range
bnmd_rel_error_total <- d |>
  filter(algorithm == "3B_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(x = iterations, y = rel_errors, color = sparsity_level, group = identifier)) +
  geom_line(alpha = .2) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(
    # subtitle = "Unequal distribution of runs & higher sparsity runs completely missing for m>=2000",
    x = "Iterations",
    y = TeX("Relative error: $\\|max\\left(0, \\Theta \\right) - X\\|_{F} / \\|X\\|_{F}$"),
    color = "Sparsity"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . , name = "Latent rank", breaks = NULL, labels = NULL)) +
  scale_x_continuous(limits = c(0, 3000), sec.axis = sec_axis(~ . , name = "Dimension", breaks = NULL, labels = NULL))

rel_error_total <- plot_grid(
  anmd_rel_error_total + theme(legend.position="none"), 
  bnmd_rel_error_total + theme(legend.position="none"), 
  labels = "AUTO")

rel_error_total_legend <- get_legend(
  # create some space to the left of the legend
  anmd_rel_error_total + theme(legend.box.margin = margin(0, 0, 0, 2))
)


rel_error_total_plots <- plot_grid(rel_error_total, rel_error_total_legend, rel_widths = c(2, .4))
rel_error_total_plots

# ANMD: first 500 iters
d |>
  filter(algorithm == "A_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(x = iterations, y = rel_errors, color = sparsity_level, group = identifier)) +
  geom_line(alpha = .2) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  # ylim(c(0,  1.25)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(
    title = "A-NMD: First 500 iterations",
    # subtitle = "Unequal distribution of runs & higher sparsity runs completely missing for m>=2000",
    x = "Iterations",
    y = TeX("Relative error: $\\|max\\left(0, \\Theta \\right) - X\\|_{F} / \\|X\\|_{F}$"),
    color = "Sparsity"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . , name = "Latent rank", breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Dimension", breaks = NULL, labels = NULL),
                     limits = c(0, 500))

# ANMD first 50 iters
d |>
  filter(algorithm == "A_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(x = iterations, y = rel_errors, color = sparsity_level, group = identifier)) +
  geom_line(alpha = .2) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(
    title = "A-NMD: First 50 iterations",
    subtitle = TeX("\\Theta0 intialized with nuclear norm"),
    x = "Iterations",
    y = TeX("Relative error: $\\|max\\left(0, \\Theta \\right) - X\\|_{F} / \\|X\\|_{F}$"),
    color = "Sparsity"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . , name = "Latent rank", breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Dimension", breaks = NULL, labels = NULL),
                     limits = c(0, 50))



# 3B NMD first 500 iters
d |>
  filter(algorithm == "3B_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(x = iterations, y = rel_errors, color = sparsity_level, group = identifier)) +
  geom_line(alpha = .2) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(
    title = "3B-NMD: First 500 iterations",
    x = "Iterations",
    y = TeX("Relative error: $\\|max\\left(0, \\Theta \\right) - X\\|_{F} / \\|X\\|_{F}$"),
    color = "Sparsity"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . , name = "Latent rank", breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Dimension", breaks = NULL, labels = NULL),
                     limits = c(0, 600))


# 3B NMD first 50 iters
d |>
  filter(algorithm == "3B_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(x = iterations, y = rel_errors, color = sparsity_level, group = identifier)) +
  geom_line(alpha = .2) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  #ylim(c(0, 2.0)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(
    title = "3B-NMD: First 50 iterations",
    subtitle = "Missing all runs for m=2000, r=8, sparsity=0.99",
    x = "Iterations",
    y = TeX("Relative error: $\\|max\\left(0, \\Theta \\right) - X\\|_{F} / \\|X\\|_{F}$"),
    color = "Sparsity"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . , name = "Latent rank", breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Dimension", breaks = NULL, labels = NULL),
                     limits = c(0, 50))
