library(ggplot2)
library(dplyr)
library(latex2exp)

folder_path <- "./simulation_nbs/lrz_square_matrices_sim_tune3B/"
csv_path <- paste0(folder_path, "all_sims_lrz_tune3B.csv")

d <- readr::read_csv(csv_path)[, -1] |>
  mutate(
    identifier = paste0(
      algorithm, m, r, sparsity_level, replication, initialization, beta
    )
  )

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

#d <- filter_by_epsilon_convergence(d)

#d <- d[d$rel_errors >= 1e-2, ]

# drop observations that hit max iteration limit
#d <- d |>
#  filter(iterations < 10000)

# d <- d |> filter(m < 5000)

# filter all
# 
rel_error_threshold <- 0.03
d <- d |>
  group_by(identifier) |> 
  mutate(converged_by_error = cumsum(cumsum(rel_errors < rel_error_threshold))) |> 
  filter(converged_by_error < 2) |> 
  ungroup() |> 
  mutate(algorithm_factor = factor(algorithm, c("A_NMD", "3B_NMD"), labels = c("A-NMD", "3B-NMD")))

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
    beta_param = unique(beta),
    replication = unique(replication),
    initialization = unique(initialization),
    final_rel_error = last(rel_errors),
    total_time = sum(times),
    mean_time_per_iter = mean(times),
    median_time_per_iter = median(times),
    n_iters = max(iterations),
    #converged_by_rel_error = max(converged_by_error)
  ) |> 
  ungroup()

# anomaly
bnmd_logtime_till_conv <- d_sum |>
filter(algorithm == "3B_NMD") |>
  mutate(sparsity_level = factor(sparsity_level),
         beta_param = factor(beta_param)) |>
  ggplot(aes(
    x = sparsity_level, y = total_time,
    color = beta_param, group = identifier
  )) + geom_jitter(alpha = .5, height = 0) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(
    x = "Sparsity",
    y = "Total time (in secs)",
    color = TeX("$\\beta$")
  )

bnmd_total_logtime_anomaly_with_grid_annot <- gridExtra::grid.arrange(
  bnmd_logtime_till_conv, top = "Dimension", right = "Rank")

ggsave("3bnmd_anomaly_multiplot.pdf", 
       plot=bnmd_total_logtime_anomaly_with_grid_annot, device = "pdf", width = 7, height = 5, 
       units = "in", dpi = 300, path = "./simulation_nbs/plots/")

# 3B whole range
d |>
  filter(algorithm == "3B_NMD", r == 32) |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(x = iterations, y = rel_errors, color = beta, group = identifier)) +
  geom_line(alpha = .2) +
  facet_grid(rows = vars(sparsity_level), cols = vars(m)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(
    title = "Preliminary: momentum 3-block NMD simulation runs",
    # subtitle = "Unequal distribution of runs & higher sparsity runs completely missing for m>=2000",
    x = "Iterations",
    y = TeX("Relative error: $\\|max\\left(0, \\Theta \\right) - X\\|_{F} / \\|X\\|_{F}$"),
    color = "beta"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . , name = "Latent rank", breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Dimension", breaks = NULL, labels = NULL))

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

# Boxplot: log(total time) over sparsity, for col: dimensions, rows: rank
d |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  mutate(r = factor(r)) |>
  group_by(identifier, algorithm, sparsity_level, m, r) |>
  summarise(total_time = sum(times)) |>
  ggplot(aes(x = r, y = log(total_time), color=algorithm)) +
  facet_grid(rows = vars(sparsity_level), cols = vars(m)) +
  geom_jitter(alpha=.5, height = 0) +
  labs(subtitle = paste0(
    "Convergence criteria:",
    "relative error < 1e-4, epsilon over 10 iterations < 1e-5,",
    "< 10 000 iterations"
  ),
  x = "Rank", y = "log(total time in secs)", fill="Algorithm")


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



# anzahl iteraationen über alle replikate
d_sum |>
  filter(algorithm == "A_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = sparsity_level, y = n_iters,
    color = sparsity_level, group = identifier
  )) + geom_jitter(alpha = .5, height = 0) +
  facet_grid(rows = vars(r), cols = vars(m))

d_sum |>
  filter(algorithm == "3B_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = sparsity_level, y = n_iters,
    color = sparsity_level, group = identifier
  )) + geom_jitter(alpha = .5, height = 0) +
  facet_grid(rows = vars(r), cols = vars(m))

d_sum |>
  filter(algorithm == "A_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(
    x = sparsity_level, y = delta_rel_error_mean,
    color = sparsity_level, group = identifier
  )) + geom_jitter(alpha = .5, height = 0) +
  facet_grid(rows = vars(r), cols = vars(m))

# ANMD whole range
d |>
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
    title = "Preliminary: Accelerated momentum NMD simulation runs",
    # subtitle = "Unequal distribution of runs & higher sparsity runs completely missing for m>=2000",
    x = "Iterations",
    y = TeX("Relative error: $\\|max\\left(0, \\Theta \\right) - X\\|_{F} / \\|X\\|_{F}$"),
    color = "Sparsity"
  )

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


# 3B whole range
d |>
  filter(algorithm == "3B_NMD") |>
  mutate(sparsity_level = factor(sparsity_level)) |>
  ggplot(aes(x = iterations, y = rel_errors, color = sparsity_level, group = identifier)) +
  geom_line(alpha = .2) +
  facet_grid(rows = vars(r), cols = vars(m)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(
    title = "Preliminary: momentum 3-block NMD simulation runs",
    # subtitle = "Unequal distribution of runs & higher sparsity runs completely missing for m>=2000",
    x = "Iterations",
    y = TeX("Relative error: $\\|max\\left(0, \\Theta \\right) - X\\|_{F} / \\|X\\|_{F}$"),
    color = "Sparsity"
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . , name = "Latent rank", breaks = NULL, labels = NULL)) +
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Dimension", breaks = NULL, labels = NULL))

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
    subtitle = "Missing all runs for m=2000, r=8, sparsity=0.99",
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
