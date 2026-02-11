# Experiment N: US-gated dVTA + No Reset
# Combines the two most promising fixes

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

# Reuse the USGated function from Exp G but with HasITI=TRUE (no reset)
# The USGated function is defined in experiment_G, let me just source that part

set.seed(42)

# Actually, let me define a minimal version inline for clarity
# (Can't source G because it runs the full experiment)

# For this test, we use the ORIGINAL Simulate.DBP but with HasITI=TRUE (no reset)
# and compare against: original+reset, original+noreset
template <- get_template_data("extinction")
n_networks <- 10

cat("=== EXPERIMENT N: Comparing Reset Strategies ===\n\n")

# Condition 1: Original with reset (baseline)
cat("--- Condition 1: Original (HasITI=FALSE, reset) ---\n")
res1 <- data.frame(M1D_t = numeric(), M1D_e = numeric())
for (net in 1:n_networks) {
  TS <- Create.Phases(phases = template$contingencies, trials = template$trials)
  r <- Simulate.DBP(NPEs = template$npes, Connections = template$connections,
    TimeSteps = TS, HasITI = c(FALSE, FALSE), disc = 0.001)
  tr <- which(r$Phase == "training"); ex <- which(r$Phase == "extinction")
  res1 <- rbind(res1, data.frame(M1D_t = r[max(tr), "M..1-D"], M1D_e = r[max(ex), "M..1-D"]))
}
cat(sprintf("  Training: %.4f  Extinction: %.4f  Change: %.1f%%\n\n",
  mean(res1$M1D_t), mean(res1$M1D_e), mean((res1$M1D_e - res1$M1D_t)/res1$M1D_t*100)))

# Condition 2: No reset for extinction only (HasITI=FALSE for training, TRUE for extinction)
cat("--- Condition 2: No reset during extinction only ---\n")
res2 <- data.frame(M1D_t = numeric(), M1D_e = numeric())
for (net in 1:n_networks) {
  TS <- Create.Phases(phases = template$contingencies, trials = template$trials)
  r <- Simulate.DBP(NPEs = template$npes, Connections = template$connections,
    TimeSteps = TS, HasITI = c(FALSE, TRUE), disc = 0.001)
  tr <- which(r$Phase == "training"); ex <- which(r$Phase == "extinction")
  res2 <- rbind(res2, data.frame(M1D_t = r[max(tr), "M..1-D"], M1D_e = r[max(ex), "M..1-D"]))
}
cat(sprintf("  Training: %.4f  Extinction: %.4f  Change: %.1f%%\n\n",
  mean(res2$M1D_t), mean(res2$M1D_e), mean((res2$M1D_e - res2$M1D_t)/res2$M1D_t*100)))

# Condition 3: No reset for extinction + higher disc
for (disc_val in c(0.005, 0.01, 0.02)) {
  cat(sprintf("--- Condition 3: No reset (ext) + disc=%.3f ---\n", disc_val))
  res3 <- data.frame(M1D_t = numeric(), M1D_e = numeric())
  for (net in 1:n_networks) {
    TS <- Create.Phases(phases = template$contingencies, trials = template$trials)
    r <- Simulate.DBP(NPEs = template$npes, Connections = template$connections,
      TimeSteps = TS, HasITI = c(FALSE, TRUE), disc = disc_val)
    tr <- which(r$Phase == "training"); ex <- which(r$Phase == "extinction")
    res3 <- rbind(res3, data.frame(M1D_t = r[max(tr), "M..1-D"], M1D_e = r[max(ex), "M..1-D"]))
  }
  pct <- mean((res3$M1D_e - res3$M1D_t)/res3$M1D_t*100)
  cat(sprintf("  Training: %.4f  Extinction: %.4f  Change: %.1f%%  >50%%ext: %d/%d\n\n",
    mean(res3$M1D_t), mean(res3$M1D_e), pct,
    sum((res3$M1D_e - res3$M1D_t)/res3$M1D_t*100 < -50), n_networks))
}

# Condition 4: No reset + ITI of 10 timesteps during extinction + disc sweep
for (disc_val in c(0.001, 0.005, 0.01)) {
  cat(sprintf("--- Condition 4: No reset + 10ts ITI + disc=%.3f ---\n", disc_val))

  res4 <- data.frame(M1D_t = numeric(), M1D_e = numeric())
  for (net in 1:n_networks) {
    trials <- template$trials
    trials$ITI <- c("US,0.00,S1,0.00,True")

    contingencies <- c(
      "training, Random, Training, 100, False",
      "extinction, Random, Extinction, 100, True, 10, 10, ITI"
    )

    TS <- Create.Phases(phases = contingencies, trials = trials)
    r <- Simulate.DBP(NPEs = template$npes, Connections = template$connections,
      TimeSteps = TS, HasITI = c(FALSE, TRUE), disc = disc_val)
    tr <- which(r$Phase == "training"); ex <- which(r$Phase == "extinction")
    res4 <- rbind(res4, data.frame(M1D_t = r[max(tr), "M..1-D"], M1D_e = r[max(ex), "M..1-D"]))
  }
  pct <- mean((res4$M1D_e - res4$M1D_t)/res4$M1D_t*100)
  cat(sprintf("  Training: %.4f  Extinction: %.4f  Change: %.1f%%  >50%%ext: %d/%d\n\n",
    mean(res4$M1D_t), mean(res4$M1D_e), pct,
    sum((res4$M1D_e - res4$M1D_t)/res4$M1D_t*100 < -50), n_networks))
}

# Show per-trial trajectory for best condition
cat("=== Best condition trajectory ===\n")
set.seed(42)
TS <- Create.Phases(phases = template$contingencies, trials = template$trials)
r <- Simulate.DBP(NPEs = template$npes, Connections = template$connections,
  TimeSteps = TS, HasITI = c(FALSE, TRUE), disc = 0.005)
tr <- which(r$Phase == "training"); ex <- which(r$Phase == "extinction")

cat("Phase     Trial  M''1->D   S1-S''1    S''1-M''1   M''1-M.1\n")
for (trial in c(1, 10, 50, 100)) {
  rows <- tr[r[tr, "Trial"] == trial]
  if (length(rows) > 0) {
    rr <- max(rows)
    cat(sprintf("Training  %3d    %.4f    %.4f     %.4f      %.4f\n",
      trial, r[rr, "M..1-D"], r[rr, "S1-S..1"], r[rr, "S..1-M..1"], r[rr, "M..1-M.1"]))
  }
}
cat("---\n")
for (trial in c(1, 5, 10, 20, 50, 100)) {
  rows <- ex[r[ex, "Trial"] == trial]
  if (length(rows) > 0) {
    rr <- max(rows)
    cat(sprintf("Extinct   %3d    %.4f    %.4f     %.4f      %.4f\n",
      trial, r[rr, "M..1-D"], r[rr, "S1-S..1"], r[rr, "S..1-M..1"], r[rr, "M..1-M.1"]))
  }
}
