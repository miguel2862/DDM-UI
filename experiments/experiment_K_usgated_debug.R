# Experiment K: Debug US-gated model - trace per-timestep weight changes

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

set.seed(42)

# We need to see the full per-timestep data including M''1 and D activations
template <- get_template_data("extinction")

# Run with ORIGINAL simulation but capture full data
TimeSteps <- Create.Phases(phases = template$contingencies, trials = template$trials)

# Use original simulation with saveData that captures all timesteps
result <- Simulate.DBP(
  NPEs = template$npes, Connections = template$connections,
  TimeSteps = TimeSteps, HasITI = template$hasITI, disc = 0.001
)

training_rows <- which(result$Phase == "training")
extinction_rows <- which(result$Phase == "extinction")

# Key columns
m1d_col <- "M..1-D"
m1_col <- "M..1"
d_col <- "D"

cat("=== Per-TS analysis: TRAINING trial 50 (CS+US) ===\n")
cat("TS  M''1->D    M''1       D          dVTA       Rule\n")
rows <- training_rows[result[training_rows, "Trial"] == 50]
for (r in rows) {
  ts <- result[r, "TimeStep"]
  w <- result[r, m1d_col]
  m1 <- result[r, m1_col]
  d <- result[r, d_col]
  dvta <- result[r, "dVTA"]
  rule <- ifelse(dvta >= 0.001, "INC", "DEC")
  cat(sprintf("%d   %.6f  %.6f  %.6f  %+.6f  %s\n", ts, w, m1, d, dvta, rule))
}

cat("\n=== Per-TS analysis: EXTINCTION trial 1 (CS only) ===\n")
rows <- extinction_rows[result[extinction_rows, "Trial"] == 1]
cat("TS  M''1->D    M''1       D          dVTA       Rule\n")
for (r in rows) {
  ts <- result[r, "TimeStep"]
  w <- result[r, m1d_col]
  m1 <- result[r, m1_col]
  d <- result[r, d_col]
  dvta <- result[r, "dVTA"]
  rule <- ifelse(dvta >= 0.001, "INC", "DEC")
  cat(sprintf("%d   %.6f  %.6f  %.6f  %+.6f  %s\n", ts, w, m1, d, dvta, rule))
}

cat("\n=== Now let's compute what the decrement SHOULD be ===\n")
cat("For M''1->D connection during extinction trial 1:\n")
cat("DECREMENT formula: Δw = -β * w * a_pre(M''1) * a_post(D)\n")
cat("β = 0.12\n\n")

beta <- 0.12
for (r in rows) {
  ts <- result[r, "TimeStep"]
  w <- result[r, m1d_col]
  m1 <- result[r, m1_col]
  d <- result[r, d_col]
  dvta <- result[r, "dVTA"]
  if (dvta < 0.001) {
    dw <- -beta * w * m1 * d
    cat(sprintf("TS%d: Δw = -0.12 * %.4f * %.4f * %.4f = %.6f (w after: %.4f)\n",
      ts, w, m1, d, dw, w + dw))
  } else {
    cat(sprintf("TS%d: INCREMENT fires (dVTA=%.4f)\n", ts, dvta))
  }
}

# Now count how many timesteps in extinction have INCREMENT vs DECREMENT
cat("\n=== INCREMENT vs DECREMENT across ALL extinction timesteps ===\n")
n_inc <- sum(result[extinction_rows, "dVTA"] >= 0.001, na.rm = TRUE)
n_dec <- sum(result[extinction_rows, "dVTA"] < 0.001, na.rm = TRUE)
total <- length(extinction_rows)
cat(sprintf("Total extinction timesteps: %d\n", total))
cat(sprintf("INCREMENT timesteps: %d (%.1f%%)\n", n_inc, n_inc/total*100))
cat(sprintf("DECREMENT timesteps: %d (%.1f%%)\n", n_dec, n_dec/total*100))

# The problem: on DECREMENT timesteps, a_post(D) is often near zero (TS1 after reset)
# Let's check D activation on decrement timesteps
dec_rows <- extinction_rows[result[extinction_rows, "dVTA"] < 0.001]
cat(sprintf("\nD activation on DECREMENT timesteps:\n"))
d_on_dec <- as.numeric(result[dec_rows, d_col])
cat(sprintf("  Mean D: %.6f\n", mean(d_on_dec, na.rm=TRUE)))
cat(sprintf("  Min D:  %.6f\n", min(d_on_dec, na.rm=TRUE)))
cat(sprintf("  Max D:  %.6f\n", max(d_on_dec, na.rm=TRUE)))
cat(sprintf("  D < 0.01: %d/%d (%.1f%%)\n",
  sum(d_on_dec < 0.01, na.rm=TRUE), length(d_on_dec),
  sum(d_on_dec < 0.01, na.rm=TRUE)/length(d_on_dec)*100))

inc_rows <- extinction_rows[result[extinction_rows, "dVTA"] >= 0.001]
cat(sprintf("\nD activation on INCREMENT timesteps:\n"))
d_on_inc <- as.numeric(result[inc_rows, d_col])
cat(sprintf("  Mean D: %.6f\n", mean(d_on_inc, na.rm=TRUE)))

cat(sprintf("\n=== THE KEY INSIGHT ===\n"))
cat("DECREMENT fires mostly on TS1 where D=0.006 (after reset)\n")
cat("Δw = -β * w * M''1 * D = -0.12 * 0.8 * 0.8 * 0.006 = -0.00046\n")
cat("INCREMENT fires on TS2-5 where D=0.7-0.9\n")
cat("The problem isn't just about when decrement fires,\n")
cat("it's that D is near-zero when decrement fires (due to reset).\n")
