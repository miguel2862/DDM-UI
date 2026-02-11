# Experiment I: No Reset + disc Sweep (DTD Adjusted approach)
#
# Key insight from experiments:
# - Exp B showed: No reset (HasITI=TRUE) gives -24.8% extinction (vs -12.7% default)
# - Exp B showed: No reset reduces %dVTA>disc from 50.7% to 10.1%
# - Exp C showed: disc=0.01 gives -23.2% with default reset
#
# Combined: No reset + disc=0.01-0.03 should give strong extinction
# while still allowing acquisition (since training has real US-driven dVTA)

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

set.seed(42)

template <- get_template_data("extinction")
n_networks <- 10

cat("=== EXPERIMENT I: No Reset + disc Sweep ===\n\n")

for (disc_val in c(0.001, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04)) {
  cat(sprintf("--- disc=%.3f, HasITI=TRUE (no reset) ---\n", disc_val))

  results <- data.frame(
    M1D_end_training = numeric(),
    M1D_end_extinction = numeric(),
    Pct_change = numeric(),
    stringsAsFactors = FALSE
  )

  for (net in 1:n_networks) {
    TimeSteps <- Create.Phases(
      phases = template$contingencies,
      trials = template$trials
    )

    # Force no reset for both phases by setting HasITI=TRUE
    result <- Simulate.DBP(
      NPEs = template$npes,
      Connections = template$connections,
      TimeSteps = TimeSteps,
      HasITI = c(TRUE, TRUE),  # No reset for either phase
      disc = disc_val
    )

    training_rows <- which(result$Phase == "training")
    extinction_rows <- which(result$Phase == "extinction")

    m1d_col <- "M..1-D"
    m1d_t <- result[max(training_rows), m1d_col]
    m1d_e <- result[max(extinction_rows), m1d_col]
    pct <- (m1d_e - m1d_t) / m1d_t * 100

    results <- rbind(results, data.frame(
      M1D_end_training = m1d_t, M1D_end_extinction = m1d_e,
      Pct_change = pct))
  }

  cat(sprintf("  Training: %.4f (SD=%.4f)  Extinction: %.4f (SD=%.4f)  Change: %.1f%% (SD=%.1f%%)\n",
    mean(results$M1D_end_training), sd(results$M1D_end_training),
    mean(results$M1D_end_extinction), sd(results$M1D_end_extinction),
    mean(results$Pct_change), sd(results$Pct_change)))
  cat(sprintf("  Nets >50%% ext: %d/%d  Nets >20%% ext: %d/%d  Nets >80%% ext: %d/%d\n\n",
    sum(results$Pct_change < -50), n_networks,
    sum(results$Pct_change < -20), n_networks,
    sum(results$Pct_change < -80), n_networks))
}

# Show per-trial trajectory for best condition
cat("\n=== Per-trial trajectory: disc=0.01, no reset, net 1 ===\n")
set.seed(42)
TimeSteps <- Create.Phases(phases = template$contingencies, trials = template$trials)
result <- Simulate.DBP(
  NPEs = template$npes, Connections = template$connections,
  TimeSteps = TimeSteps, HasITI = c(TRUE, TRUE), disc = 0.01
)

training_rows <- which(result$Phase == "training")
extinction_rows <- which(result$Phase == "extinction")
m1d_col <- "M..1-D"

cat("Phase     Trial  M''1->D\n")
for (trial in c(1, 10, 20, 50, 80, 100)) {
  rows <- training_rows[result[training_rows, "Trial"] == trial]
  if (length(rows) > 0) cat(sprintf("Training  %3d    %.4f\n", trial, result[max(rows), m1d_col]))
}
cat("---\n")
for (trial in c(1, 5, 10, 20, 30, 50, 70, 100)) {
  rows <- extinction_rows[result[extinction_rows, "Trial"] == trial]
  if (length(rows) > 0) cat(sprintf("Extinct   %3d    %.4f\n", trial, result[max(rows), m1d_col]))
}

# Show all connection weights
cat("\n--- All connection weights ---\n")
conn_cols <- grep("-", colnames(result), value = TRUE)
cat(sprintf("%-20s %12s %12s %12s\n", "Connection", "End_train", "End_ext", "Change"))
for (col in conn_cols) {
  train_val <- result[max(training_rows), col]
  ext_val <- result[max(extinction_rows), col]
  cat(sprintf("%-20s %12.6f %12.6f %12.6f\n", col, train_val, ext_val, ext_val - train_val))
}
