# Experiment L: Quick trace of US-gated extinction

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

set.seed(42)

template <- get_template_data("extinction")
TimeSteps <- Create.Phases(phases = template$contingencies, trials = template$trials)

# Run original
result_orig <- Simulate.DBP(
  NPEs = template$npes, Connections = template$connections,
  TimeSteps = TimeSteps, HasITI = template$hasITI, disc = 0.001
)

extinction_rows <- which(result_orig$Phase == "extinction")
m1d_col <- "M..1-D"
d_col <- "D"

# Count: what % of decrement timesteps have D > 0.5?
dec_rows <- extinction_rows[result_orig[extinction_rows, "dVTA"] < 0.001]
d_on_dec <- as.numeric(result_orig[dec_rows, d_col])
cat("=== ORIGINAL: D activation when DECREMENT fires ===\n")
cat(sprintf("Total decrement timesteps: %d\n", length(dec_rows)))
cat(sprintf("D > 0.5: %d (%.1f%%)\n", sum(d_on_dec > 0.5, na.rm=TRUE), sum(d_on_dec > 0.5, na.rm=TRUE)/length(d_on_dec)*100))
cat(sprintf("D > 0.1: %d (%.1f%%)\n", sum(d_on_dec > 0.1, na.rm=TRUE), sum(d_on_dec > 0.1, na.rm=TRUE)/length(d_on_dec)*100))
cat(sprintf("D < 0.01: %d (%.1f%%)\n", sum(d_on_dec < 0.01, na.rm=TRUE), sum(d_on_dec < 0.01, na.rm=TRUE)/length(d_on_dec)*100))

# For each decrement timestep, compute the actual decrement
beta <- 0.12
cat("\nExpected decrement magnitudes on non-reset timesteps:\n")
cat("TS  D        M''1     w        Δw\n")
for (r in dec_rows[1:min(10, length(dec_rows))]) {
  d <- as.numeric(result_orig[r, d_col])
  m1 <- as.numeric(result_orig[r, "M..1"])
  w <- as.numeric(result_orig[r, m1d_col])
  ts <- result_orig[r, "TimeStep"]
  trial <- result_orig[r, "Trial"]
  dw <- -beta * w * m1 * d
  cat(sprintf("T%d.%d  %.4f  %.4f  %.4f  %+.6f\n", trial, ts, d, m1, w, dw))
}

# KEY: the issue with US-gated is that dCA1 also depends on dVTA
# When dVTA=0, dCA1 = dH + 0*(1-prev) = just dH
# For sensory/hippocampal connections, d = dCA1
# If dCA1 is still positive (due to dH component), those connections get INCREMENT
cat("\n=== dCA1 ANALYSIS ===\n")
cat("In US-gated model, dVTA=0 during extinction.\n")
cat("dCA1 = dH + dVTA*(1-prev_dCA1) = dH + 0 = dH\n")
cat("dH = mean(|H(t) - H(t-1)|) for hippocampal units\n\n")

# Check dH values during extinction
dh_vals <- as.numeric(result_orig[extinction_rows, "dH"])
cat(sprintf("dH during extinction: mean=%.6f, min=%.6f, max=%.6f\n",
  mean(dh_vals, na.rm=TRUE), min(dh_vals, na.rm=TRUE), max(dh_vals, na.rm=TRUE)))
cat(sprintf("dH > 0.001 (i.e., increment for sensory layers): %d/%d (%.1f%%)\n",
  sum(dh_vals > 0.001, na.rm=TRUE), length(dh_vals), sum(dh_vals > 0.001, na.rm=TRUE)/length(dh_vals)*100))

# The issue: S1->S''1 connection uses dCA1
# If dCA1 > disc, S1->S''1 gets INCREMENT
# This keeps S''1 active, which keeps M''1 active, which keeps D active...
# And for M''1->D, which uses dVTA (dopaminergic layer), dVTA=0 in US-gated
# So M''1->D should get pure decrement

# BUT WAIT - let me check what a_pre means for M''1->D
# The connection is M''1 -> D
# pre = M''1, post = D
# In the decrement formula: Δw = -β * w * a_pre * a_post
# a_pre = activation of M''1 (presynaptic)
# a_post = activation of D (postsynaptic)
# But in the code, "network[[i]]" is the postsynaptic unit (D)
# and network[[pre]] is the presynaptic unit (M''1)

# Let me trace exactly what happens to M''1->D in the first 5 extinction trials
# by examining the actual weight at each timestep
cat("\n=== M''1->D weight per timestep in extinction ===\n")
cat("Trial  TS  Weight\n")
for (trial in 1:5) {
  rows <- extinction_rows[result_orig[extinction_rows, "Trial"] == trial]
  for (r in rows) {
    ts <- result_orig[r, "TimeStep"]
    w <- result_orig[r, m1d_col]
    cat(sprintf("%3d    %d   %.6f\n", trial, ts, w))
  }
  cat("---\n")
}
