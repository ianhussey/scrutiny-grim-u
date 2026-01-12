library(tidyverse)
library(roundwork)

# core engine
grimu_map_pvalues <- function(n1, n2, u_min = NULL, u_max = NULL, alternative = "two.sided") {
  
  # Safety Check: Input Validation
  if (is.na(n1) || is.na(n2)) return(tibble(U = numeric(), is_integer = logical()))
  
  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  N <- n1 + n2
  mu <- (n1 * n2) / 2
  max_u <- n1 * n2
  
  # Default bounds logic (Dynamic Tail Selection)
  if (is.null(u_min) || is.null(u_max)) {
    if (alternative == "greater") {
      if (is.null(u_min)) u_min <- floor(mu) 
      if (is.null(u_max)) u_max <- max_u
    } else {
      if (is.null(u_min)) u_min <- 0
      if (is.null(u_max)) u_max <- ceiling(mu) 
    }
  }
  
  # Enforce Half-Integer Lattice
  # Snap arbitrary bounds to the nearest valid U values (integers or half-integers).
  # If u_min = 0.2 -> snaps up to 0.5
  # If u_max = 0.8 -> snaps down to 0.5
  
  # Handle NULLs before math (though logic above ensures they aren't NULL)
  if (is.null(u_min)) u_min <- 0
  if (is.null(u_max)) u_max <- max_u
  
  u_start <- ceiling(u_min * 2) / 2
  u_end   <- floor(u_max * 2) / 2
  
  # Safety Clamp: Ensure bounds are physical
  u_start <- max(0, u_start)
  u_end   <- min(max_u, u_end)
  
  # Safety Check: If start > end, return empty tibble immediately
  if (u_start > u_end) return(tibble(U = numeric(), is_integer = logical()))
  
  vals <- roundwork::round_up(seq(u_start, u_end, by = 0.5), 1)
  
  # Constants
  sigma_no_ties <- sqrt((n1 * n2 * (N + 1)) / 12)
  correction_term <- (n1 * n2 * 6) / (12 * N * (N - 1))
  sigma_one_tie <- sqrt((n1 * n2 * (N + 1)) / 12 - correction_term)
  
  # Vectorized Calc
  results_df <- tibble(U = vals) %>%
    mutate(
      is_integer = (U %% 1 == 0),
      
      # Exact
      p_exact = if_else(is_integer, case_when(
        alternative == "less"      ~ pwilcox(U, n1, n2),
        alternative == "greater"   ~ pwilcox(U - 1, n1, n2, lower.tail = FALSE),
        alternative == "two.sided" ~ {
          p_lower <- pwilcox(U, n1, n2)
          p_upper <- pwilcox(U - 1, n1, n2, lower.tail = FALSE)
          2 * pmin(p_lower, p_upper)
        }
      ), NA_real_),
      
      # Deviations
      dev_cc = case_when(
        alternative == "two.sided" ~ pmax(0, abs(U - mu) - 0.5),
        alternative == "less"      ~ (U - mu) + 0.5,
        alternative == "greater"   ~ (U - mu) - 0.5
      ),
      dev_uncorr = case_when(
        alternative == "two.sided" ~ abs(U - mu),
        alternative == "less"      ~ (U - mu),
        alternative == "greater"   ~ (U - mu)
      ),
      
      # P-values
      # No Ties
      z_corr_no_ties = dev_cc / sigma_no_ties,
      p_corr_no_ties = case_when(
        alternative == "two.sided" ~ 2 * pnorm(z_corr_no_ties, lower.tail = FALSE),
        alternative == "less"      ~ pnorm(z_corr_no_ties, lower.tail = TRUE),
        alternative == "greater"   ~ pnorm(z_corr_no_ties, lower.tail = FALSE)
      ),
      
      # No Ties (Uncorrected)
      z_uncorr_no_ties = dev_uncorr / sigma_no_ties,
      p_uncorr_no_ties = case_when(
        alternative == "two.sided" ~ 2 * pnorm(z_uncorr_no_ties, lower.tail = FALSE),
        alternative == "less"      ~ pnorm(z_uncorr_no_ties, lower.tail = TRUE),
        alternative == "greater"   ~ pnorm(z_uncorr_no_ties, lower.tail = FALSE)
      ),
      
      # Tied (Corrected)
      z_corr_tied = dev_cc / sigma_one_tie,
      p_corr_tied = case_when(
        alternative == "two.sided" ~ 2 * pnorm(z_corr_tied, lower.tail = FALSE),
        alternative == "less"      ~ pnorm(z_corr_tied, lower.tail = TRUE),
        alternative == "greater"   ~ pnorm(z_corr_tied, lower.tail = FALSE)
      ),
      
      # Tied (Uncorrected)
      z_uncorr_tied = dev_uncorr / sigma_one_tie,
      p_uncorr_tied = case_when(
        alternative == "two.sided" ~ 2 * pnorm(z_uncorr_tied, lower.tail = FALSE),
        alternative == "less"      ~ pnorm(z_uncorr_tied, lower.tail = TRUE),
        alternative == "greater"   ~ pnorm(z_uncorr_tied, lower.tail = FALSE)
      )
    ) %>%
    # Final Clamp: Ensure p <= 1 (Standard behavior)
    mutate(across(starts_with("p_"), ~ pmin(1, .))) %>%
    select(U, is_integer, starts_with("p_"))
  
  return(results_df)
}

# forensic tool
grimu_check <- function(n1, n2, p_reported, 
                        comparison = "equal", 
                        digits = 2, 
                        rounding = NULL, 
                        p_min = NULL, p_max = NULL, 
                        alternative = "two.sided") {
  
  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  
  # --- 1. Range Detection ---
  if (is.null(p_min) || is.null(p_max)) {
    # Default window for bounds search (wide enough to catch all rounding types)
    window <- 10^(-digits) * 5 
    p_max_search <- p_reported + window
    
    # Calculate raw lower bound
    raw_p_min <- p_reported - window
    
    # If comparison is inequality, or window touches 0, anchor the "deep tail".
    if (comparison == "less_than" || raw_p_min <= 0) {
      p_min_search <- NA_real_ 
    } else {
      # Standard case: Window is strictly positive (e.g., 0.04 to 0.06)
      p_min_search <- raw_p_min
    }
  } else {
    p_min_search <- p_min
    p_max_search <- p_max
    # If user manually provides 0, treat it as NA for Z-score purposes
    if (p_min_search <= 0) p_min_search <- NA_real_
  }
  
  # --- 2. U Bounds Calculation ---
  N <- n1 + n2
  mu <- (n1 * n2) / 2
  sigma_est <- sqrt((n1 * n2 * (N + 1)) / 12)
  max_u <- n1 * n2
  
  # Helper: Convert P to Z (Signed)
  p_to_z <- function(p, alt) {
    p_safe <- min(1, max(0, p)) # clamp p to [0,1]
    if (is.na(p)) return(Inf) # Infinite Z implies maximal U deviation
    if (alt == "two.sided") return(qnorm(1 - p / 2)) # Always positive distance
    if (alt == "less")      return(qnorm(p))         # Negative Z
    if (alt == "greater")   return(qnorm(1 - p))     # Positive Z
  }
  
  # Calculate Z boundaries for the p-value window
  # Note: 
  # Small P -> Large Z magnitude (Deep Tail)
  # Large P -> Small Z magnitude (Near Mean)
  
  z_deep <- p_to_z(if(is.na(p_min_search)) 0 else p_min_search, alternative)
  z_shallow <- p_to_z(p_max_search, alternative)
  
  # Convert Z to U
  # U ~ mu + Z*sigma
  u_bound_1 <- mu + z_deep * sigma_est
  u_bound_2 <- mu + z_shallow * sigma_est
  
  # NA Safety: If bounds are NaN (e.g., extremely far out), default to full range
  if (is.na(u_bound_1) || is.nan(u_bound_1)) u_bound_1 <- if(alternative=="greater") max_u else 0
  if (is.na(u_bound_2) || is.nan(u_bound_2)) u_bound_2 <- if(alternative=="greater") mu else mu
  
  # Sort bounds
  raw_start <- min(u_bound_1, u_bound_2)
  raw_end   <- max(u_bound_1, u_bound_2)
  
  # Pad and Clamp
  u_start_est <- floor(raw_start - 2)
  u_end_est   <- ceiling(raw_end + 2)
  
  # Physical Clamping
  u_start_est <- max(0, u_start_est)
  u_end_est   <- min(max_u, u_end_est)
  
  # --- 3. Call Engine ---
  results_df <- grimu_map_pvalues(n1, n2, u_min = u_start_est, u_max = u_end_est, 
                                  alternative = alternative)
  
  # --- 4. Consistency Logic (Interval Union) ---
  if (is.null(rounding)) {
    methods <- c("round", "trunc") 
  } else {
    methods <- match.arg(rounding, c("round", "trunc", "up"), several.ok = TRUE)
  }
  
  epsilon <- 10^(-digits)
  buffer <- 1e-14 
  
  is_match <- function(val) {
    if (is.na(val)) return(FALSE)
    if ("round" %in% methods) {
      if (val >= (p_reported - 0.5 * epsilon - buffer) && 
          val <  (p_reported + 0.5 * epsilon - buffer)) return(TRUE)
    }
    if ("trunc" %in% methods) {
      if (val >= (p_reported - buffer) && 
          val <  (p_reported + epsilon - buffer)) return(TRUE)
    }
    if ("up" %in% methods) {
      if (val >  (p_reported - epsilon + buffer) && 
          val <= (p_reported + buffer)) return(TRUE)
    }
    return(FALSE)
  }
  
  check_col <- function(col_val) {
    if (comparison == "equal") {
      return(is_match(col_val))
    } else {
      return(!is.na(col_val) && col_val < p_reported)
    }
  }
  
  results_checked <- results_df %>%
    rowwise() %>%
    mutate(
      valid_exact          = check_col(p_exact),
      valid_corr_no_ties   = check_col(p_corr_no_ties),
      valid_uncorr_no_ties = check_col(p_uncorr_no_ties),
      valid_corr_tied      = check_col(p_corr_tied),
      valid_uncorr_tied    = check_col(p_uncorr_tied),
      
      is_consistent = valid_exact | 
        valid_corr_no_ties | valid_uncorr_no_ties | 
        valid_corr_tied | valid_uncorr_tied
    ) %>%
    ungroup() %>%
    # Filter for display (Diagnostic mode)
    filter(
      is_consistent | 
        (p_corr_tied >= (if(is.na(p_min_search)) 0 else p_min_search) & 
           p_corr_tied <= p_max_search)
    )
  
  summary_df <- tibble(
    n1 = n1, n2 = n2, p_reported = p_reported,
    alternative = alternative,
    rounding = paste(methods, collapse = "+"),
    consistent = any(results_checked$is_consistent),
    matches_exact = any(results_checked$valid_exact),
    matches_no_ties = any(results_checked$valid_corr_no_ties | results_checked$valid_uncorr_no_ties),
    matches_ties = any(results_checked$valid_corr_tied | results_checked$valid_uncorr_tied)
  )
  
  return(list(summary = summary_df, details = results_checked))
}

# saturation calc
grimu_saturation <- function(n1, n2, decimals = 3, p_lower_threshold = 0, p_upper_threshold = 1) {
  
  # 1. CALL THE ENGINE 
  # Get all U values from 0 to Mean (Sufficient because distribution is symmetric)
  # This covers every possible unique p-value the test can produce.
  mu <- (n1 * n2) / 2
  p_space <- grimu_map_pvalues(n1, n2, u_min = 0, u_max = floor(mu))
  
  # 2. Reshape and Filter
  unique_rounded_p <- p_space %>%
    pivot_longer(cols = starts_with("p_"), values_to = "p_val") %>%
    filter(!is.na(p_val)) %>%
    filter(p_val >= p_lower_threshold & p_val <= p_upper_threshold) %>%
    # Round to target precision
    mutate(p_rounded = roundwork::round_up(p_val, decimals)) %>%
    distinct(p_rounded) %>%
    nrow()
  
  # 3. Calculate Coverage
  total_slots <- length(seq(p_lower_threshold, p_upper_threshold, by = 10^-decimals))
  
  return(unique_rounded_p / total_slots)
}
