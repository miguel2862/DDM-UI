# Experiment J: Debug why US-gated dVTA decrement is so weak
# dVTA=0 should put d<disc on EVERY timestep
# So decrement should fire on ALL 500 timesteps
# β * w * a_pre * a_post should give decent decrements

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

set.seed(42)

template <- get_template_data("extinction")
TimeSteps <- Create.Phases(phases = template$contingencies, trials = template$trials)

# Run original and check per-timestep dVTA and weight changes during extinction
result_orig <- Simulate.DBP(
  NPEs = template$npes, Connections = template$connections,
  TimeSteps = TimeSteps, HasITI = template$hasITI, disc = 0.001
)

training_rows <- which(result_orig$Phase == "training")
extinction_rows <- which(result_orig$Phase == "extinction")
m1d_col <- "M..1-D"

cat("=== Original Model: Per-Timestep Extinction Analysis ===\n\n")
cat("Trial  TS  M''1->D         dVTA        d>=disc?  Rule\n")
for (trial in c(1, 2, 50, 100)) {
  rows <- extinction_rows[result_orig[extinction_rows, "Trial"] == trial]
  for (r in rows) {
    m1d <- result_orig[r, m1d_col]
    dvta <- result_orig[r, "dVTA"]
    ts <- result_orig[r, "TimeStep"]
    rule <- ifelse(dvta >= 0.001, "INCREMENT", "DECREMENT")
    cat(sprintf("%3d    %d   %.6f      %+.6f    %s      %s\n", trial, ts, m1d, dvta, ifelse(dvta >= 0.001, "YES", "NO "), rule))
  }
  cat("---\n")
}

# Calculate what SHOULD happen with pure decrement:
# If decrement fires every timestep during extinction:
# Δw = -β * w * a_pre * a_post per timestep
# β=0.12, a_pre(M''1)≈0.8, a_post(D)≈0.8
# After 1 trial (5 ts): w(n+1) = w(n) * (1 - 0.12 * 0.8 * 0.8)^5 = w(n) * 0.923^5 = w(n) * 0.674
# After 100 trials: w * 0.674^100 ≈ essentially 0!

cat("\n=== Theoretical Exponential Decrement Analysis ===\n")
w0 <- 0.825  # typical end-of-training weight
beta <- 0.12
a_pre <- 0.8  # typical M''1 activation
a_post <- 0.8  # typical D activation
decay_per_ts <- beta * a_pre * a_post
cat(sprintf("β=%.2f, a_pre=%.2f, a_post=%.2f\n", beta, a_pre, a_post))
cat(sprintf("Decay rate per TS: β*a_pre*a_post = %.4f\n", decay_per_ts))
cat(sprintf("After 1 TS: w = %.4f * (1 - %.4f) = %.4f\n", w0, decay_per_ts, w0*(1-decay_per_ts)))
cat(sprintf("After 5 TS: w = %.4f * (1 - %.4f)^5 = %.4f\n", w0, decay_per_ts, w0*(1-decay_per_ts)^5))
cat(sprintf("After 25 TS (5 trials): w = %.4f\n", w0*(1-decay_per_ts)^25))
cat(sprintf("After 500 TS (100 trials): w = %.8f\n", w0*(1-decay_per_ts)^500))

# But the actual decrement uses current weight * β * a_pre * a_post
# which is exponential decay - let me check what the actual a_pre and a_post are
# during extinction in the US-gated model

# Now: What are the ACTUAL a_pre (presynaptic = M''1) and a_post (postsynaptic = D)
# during extinction?
cat("\n=== Actual M''1 and D activations during extinction ===\n")
cat("Trial  TS  M''1          D\n")
m1_col <- "M..1"
d_col <- "D"
for (trial in c(1, 2, 10, 50, 100)) {
  rows <- extinction_rows[result_orig[extinction_rows, "Trial"] == trial]
  for (r in rows) {
    m1_act <- result_orig[r, m1_col]
    d_act <- result_orig[r, d_col]
    ts <- result_orig[r, "TimeStep"]
    cat(sprintf("%3d    %d   %.6f      %.6f\n", trial, ts, m1_act, d_act))
  }
  cat("---\n")
}
