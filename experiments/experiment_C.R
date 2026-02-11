# Experiment C: disc criterion sweep with more granularity
# Tests disc = 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

set.seed(42)

template <- get_template_data("extinction")
n_networks <- 5
disc_values <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5)

cat("=== EXPERIMENT C: disc Criterion Sweep ===\n\n")

results_C <- data.frame(
  disc = numeric(),
  Network = integer(),
  M1D_end_training = numeric(),
  M1D_end_extinction = numeric(),
  M1D_change = numeric(),
  Pct_change = numeric(),
  Pct_dVTA_above_disc = numeric(),
  stringsAsFactors = FALSE
)

for (disc_val in disc_values) {
  cat(sprintf("--- disc = %.3f ---\n", disc_val))

  for (net in 1:n_networks) {
    TimeSteps <- Create.Phases(
      phases = template$contingencies,
      trials = template$trials
    )

    result <- Simulate.DBP(
      NPEs = template$npes,
      Connections = template$connections,
      TimeSteps = TimeSteps,
      HasITI = template$hasITI,
      disc = disc_val
    )

    training_rows <- which(result$Phase == "training")
    extinction_rows <- which(result$Phase == "extinction")

    m1d_col <- "M..1-D"
    m1d_end_training <- result[max(training_rows), m1d_col]
    m1d_end_extinction <- result[max(extinction_rows), m1d_col]
    m1d_change <- m1d_end_extinction - m1d_end_training
    pct_change <- m1d_change / m1d_end_training * 100

    ext_dvta <- result[extinction_rows, "dVTA"]
    pct_above <- sum(ext_dvta > disc_val, na.rm = TRUE) / length(ext_dvta) * 100

    results_C <- rbind(results_C, data.frame(
      disc = disc_val,
      Network = net,
      M1D_end_training = round(m1d_end_training, 6),
      M1D_end_extinction = round(m1d_end_extinction, 6),
      M1D_change = round(m1d_change, 6),
      Pct_change = round(pct_change, 2),
      Pct_dVTA_above_disc = round(pct_above, 2),
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  Net %d: M''1->D train=%.4f, ext=%.4f, change=%.4f (%.1f%%), %%dVTA>disc=%.1f%%\n",
                net, m1d_end_training, m1d_end_extinction, m1d_change, pct_change, pct_above))
  }
}

cat("\n=== EXPERIMENT C SUMMARY ===\n")
cat(sprintf("\n%-8s %-12s %-12s %-12s %-12s %-12s %-12s\n",
            "disc", "M_train", "M_ext", "Change", "Pct_chg", "%dVTA>disc", "#>50%dec"))

for (disc_val in disc_values) {
  subset <- results_C[results_C$disc == disc_val, ]
  n_gt50 <- sum(subset$Pct_change < -50)

  cat(sprintf("%-8.3f %-12.4f %-12.4f %-12.4f %-12.1f %-12.1f %-12d\n",
              disc_val,
              mean(subset$M1D_end_training),
              mean(subset$M1D_end_extinction),
              mean(subset$M1D_change),
              mean(subset$Pct_change),
              mean(subset$Pct_dVTA_above_disc),
              n_gt50))
}
