library(dplyr)

print_removed_rows <- function(before, after) {
  difference <- before - after
  print(paste0("Removed ", difference, " rows."))
}

filter_by_maxiters <- function(df, maxiters = 10000) {
  # remove all runs that reached maxiters
  n_before <- nrow(df)
  n_runs_before <- length(unique(df$identifier))
  
  result <- df |> 
    group_by(identifier) |> 
    mutate(max_iters = max(iterations)) |> 
    filter(max_iters < 10000) |> 
    select(-max_iters) |> 
    ungroup()
  
  print(paste0("Removed ", n_runs_before - length(unique(result$identifier)), " simulation runs."))
  result
}

filter_by_final_relerror <- function(df, final_rel_error = 1e-5) {
  # remove all iterations after final_rel_error was reached
  n_before <- nrow(df)
  result <- df |> 
    filter(rel_errors >= final_rel_error)
  print_removed_rows(n_before, nrow(result))
  result
}

filter_by_epsilon_over_iters <- function(df, epsilon = 1e-5, over_n_iters = 10) {
  # remove all iterations after gradient epsilon is reached over `over_n_iters`
  n_before <- nrow(df)
  result <- df |>
    group_by(m, r, sparsity_level, replication, initialization) |>
    mutate(
      has_converged_epsilon = abs(
        rel_errors - lag(rel_errors, over_n_iters)
      ) < epsilon
    ) |>
    ungroup() |>
    mutate(
      has_converged_epsilon = tidyr::replace_na(has_converged_epsilon, FALSE)
    ) |>
    filter(!has_converged_epsilon) |> 
    select(-has_converged_epsilon)
  print_removed_rows(n_before, nrow(result))
  result
}

make_summary_df <- function(df) {
  d_nonvarying <- df |> 
    select(norm_X:initialization, algorithm, identifier) |> 
    distinct()
  
  d_varying <- df |> 
    select(identifier, times, rel_errors, iterations) |> 
    group_by(identifier) |> 
    mutate(delta_rel_error_last_iter = rel_errors - lag(rel_errors)) |> 
    summarize(
      time_total = sum(times),
      n_iters = max(iterations),
      time_per_iter_mean = mean(times),
      time_per_iter_median = median(times),
      init_rel_error = first(rel_errors),
      final_rel_error = last(rel_errors),
      delta_rel_error_mean = mean(delta_rel_error_last_iter, na.rm = TRUE))
  
  d_nonvarying |> left_join(d_varying, by = "identifier")
}
