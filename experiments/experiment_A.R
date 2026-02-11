# Experiment A: Varying timesteps per trial (5, 8, 12, 20)
# Tests how the number of timesteps per extinction trial affects weight change

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

set.seed(42)

template <- get_template_data("extinction")

# Define extinction trials with different numbers of timesteps
make_extinction_trials <- function(n_ts) {
  replicate(n_ts, "US,0.00,S1,1.00,True", simplify = TRUE)
}

timestep_counts <- c(5, 8, 12, 20)
n_networks <- 5

cat("=== EXPERIMENT A: Varying Timesteps Per Extinction Trial ===\n\n")

results_A <- data.frame(
  Timesteps = integer(),
  Network = integer(),
  M1D_end_training = numeric(),
  M1D_end_extinction = numeric(),
  M1D_change = numeric(),
  Pct_dVTA_above_disc = numeric(),
  stringsAsFactors = FALSE
)

for (n_ts in timestep_counts) {
  cat(sprintf("--- Testing %d timesteps per trial ---\n", n_ts))

  for (net in 1:n_networks) {
    # Training trials stay the same (5 timesteps, US on last)
    trials <- template$trials
    # Modify extinction trial to have n_ts timesteps (all CS only)
    trials$Extinction <- make_extinction_trials(n_ts)

    # Contingencies: training stays same, extinction uses modified trials
    contingencies <- template$contingencies

    TimeSteps <- Create.Phases(
      phases = contingencies,
      trials = trials
    )

    result <- Simulate.DBP(
      NPEs = template$npes,
      Connections = template$connections,
      TimeSteps = TimeSteps,
      HasITI = template$hasITI,
      disc = 0.001
    )

    # Get M''1->D weight at end of training and end of extinction
    training_rows <- which(result$Phase == "training")
    extinction_rows <- which(result$Phase == "extinction")

    m1d_col <- "M..1-D"
    m1d_end_training <- result[max(training_rows), m1d_col]
    m1d_end_extinction <- result[max(extinction_rows), m1d_col]
    m1d_change <- m1d_end_extinction - m1d_end_training

    # % of extinction timesteps with dVTA > disc
    ext_dvta <- result[extinction_rows, "dVTA"]
    pct_above <- sum(ext_dvta > 0.001, na.rm = TRUE) / length(ext_dvta) * 100

    results_A <- rbind(results_A, data.frame(
      Timesteps = n_ts,
      Network = net,
      M1D_end_training = round(m1d_end_training, 6),
      M1D_end_extinction = round(m1d_end_extinction, 6),
      M1D_change = round(m1d_change, 6),
      Pct_dVTA_above_disc = round(pct_above, 2),
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  Net %d: M''1->D training=%.4f, extinction=%.4f, change=%.4f, %%dVTA>disc=%.1f%%\n",
                net, m1d_end_training, m1d_end_extinction, m1d_change, pct_above))
  }
}

cat("\n=== EXPERIMENT A SUMMARY ===\n")
for (n_ts in timestep_counts) {
  subset <- results_A[results_A$Timesteps == n_ts, ]
  cat(sprintf("\n%d timesteps/trial:\n", n_ts))
  cat(sprintf("  Mean M''1->D end training:    %.4f (SD=%.4f)\n", mean(subset$M1D_end_training), sd(subset$M1D_end_training)))
  cat(sprintf("  Mean M''1->D end extinction:  %.4f (SD=%.4f)\n", mean(subset$M1D_end_extinction), sd(subset$M1D_end_extinction)))
  cat(sprintf("  Mean change:                  %.4f (SD=%.4f)\n", mean(subset$M1D_change), sd(subset$M1D_change)))
  cat(sprintf("  Mean %%dVTA > disc:            %.1f%% (SD=%.1f%%)\n", mean(subset$Pct_dVTA_above_disc), sd(subset$Pct_dVTA_above_disc)))
  cat(sprintf("  Pct change from training:     %.1f%%\n", mean(subset$M1D_change) / mean(subset$M1D_end_training) * 100))
}
