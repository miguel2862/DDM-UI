# Experiment M: Full DTD Adjusted Fix
#
# Two corrections:
# 1. US-gated dVTA: dVTA = D_actual - D_conditioned (only US surplus)
# 2. US-gated dCA1: dCA1 = dH + dVTA*(1-prev) only when dVTA > 0
#    When dVTA = 0 (no US), dCA1 = 0 too (no learning signal for sensory layers)
#
# This means during extinction (no US):
# - dVTA = 0 → DECREMENT fires for motor/dopaminergic connections
# - dCA1 = 0 → DECREMENT fires for sensory/hippocampal connections
# - ALL connections get pure decrement during extinction

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

set.seed(42)

Simulate.DBP.FullFix <- function(NPEs, Connections, TimeSteps, HasITI = c(FALSE),
                                   disc = 0.001, threshold = "gaussian",
                                   pupdate = "async_random",
                                   saveData = list()) {

  tryCatch({
    setClass("NPE", slots = c(
      Activation = "numeric", PreviousActivation = "numeric",
      ActivationDecay = "numeric", ExcitatoryInput = "numeric",
      InhibitoryInput = "numeric", PreviousExcitatoryInput = "numeric",
      TemporalSummation = "numeric", Name = "character",
      Type = "character", Layer = "character", Threshold = "numeric",
      mu = "numeric", sigma = "numeric", logisSigma = "numeric",
      InputConnections = "list", r = "numeric"))
    setClass("Connection", slots = c(
      weight = "numeric", Name = "character", alpha = "numeric",
      beta = "numeric", alpha_prime = "numeric", beta_prime = "numeric",
      preSinapticNPE = "character", p = "numeric"))

    NPEs <- as.data.frame(NPEs, stringsAsFactors = FALSE)
    Connections <- as.data.frame(Connections, stringsAsFactors = FALSE)
    NPEs[, 1:3] <- lapply(NPEs[, 1:3, drop = FALSE], as.character)
    NPEs[, 4:9] <- lapply(NPEs[, 4:9, drop = FALSE], function(x) suppressWarnings(as.numeric(x)))
    Connections[, 1:2] <- lapply(Connections[, 1:2, drop = FALSE], as.character)
    Connections[, 3:7] <- lapply(Connections[, 3:7, drop = FALSE], function(x) suppressWarnings(as.numeric(x)))

    ComputeInputs <- function(npe) {
      inputs <- c(0, 0)
      if (length(npe@InputConnections) > 0) {
        for (j in 1:length(npe@InputConnections)) {
          pre <- npe@InputConnections[[j]]@preSinapticNPE
          if (network[[pre]]@Type == "Excitatory") {
            inputs[1] <- inputs[1] + network[[pre]]@Activation * npe@InputConnections[[j]]@weight
          } else {
            inputs[2] <- inputs[2] + network[[pre]]@Activation * npe@InputConnections[[j]]@weight
          }
        }
      }
      return(inputs)
    }
    L <- function(x, sigma) { return(1 / (1 + exp((-x + 0.5) / sigma))) }
    estBetaParams <- function(mu, sigma) {
      var <- sigma^2; alpha <- mu * (mu - mu^2 - var) / var
      beta <- (mu * (1 - mu) / var - 1) * (1 - mu)
      if (alpha <= 0 | beta <= 0) stop("invalid beta params")
      return(params = list(alpha = alpha, beta = beta))
    }

    # FIX 1: US-gated dVTA
    dVTA_us_gated <- function() {
      d_actual <- 0; d_conditioned <- 0; n <- 0
      for (i in 1:length(network)) {
        if (network[[i]]@Layer == "Dopaminergic") {
          d_actual <- d_actual + network[[i]]@Activation; n <- n + 1
          us_driving <- FALSE
          for (j in 1:length(network[[i]]@InputConnections)) {
            pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
            if (network[[pre]]@Layer == "US" && network[[pre]]@Activation > 0) { us_driving <- TRUE; break }
          }
          if (us_driving) {
            exc_no_us <- 0; inh <- 0
            for (j in 1:length(network[[i]]@InputConnections)) {
              pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
              if (network[[pre]]@Layer != "US") {
                if (network[[pre]]@Type == "Excitatory") {
                  exc_no_us <- exc_no_us + network[[pre]]@Activation * network[[i]]@InputConnections[[j]]@weight
                } else { inh <- inh + network[[pre]]@Activation * network[[i]]@InputConnections[[j]]@weight }
              }
            }
            p_e <- L(exc_no_us, network[[i]]@logisSigma); p_i <- L(inh, network[[i]]@logisSigma)
            if (p_e > p_i && p_e >= network[[i]]@Threshold) {
              d_conditioned <- d_conditioned + p_e + network[[i]]@TemporalSummation *
                L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) * (1 - p_e) - p_i
            } else if (p_e > p_i) {
              d_conditioned <- d_conditioned + L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) -
                network[[i]]@ActivationDecay * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma)
            }
          } else { d_conditioned <- d_conditioned + network[[i]]@Activation }
        }
      }
      if (n == 0) return(0)
      return((d_actual - d_conditioned) / n)
    }

    # FIX 2: US-gated dCA1
    # When dVTA > 0 (US present): dCA1 = dH + dVTA*(1-prev) [original]
    # When dVTA <= 0 (no US): dCA1 = 0 [forces decrement for sensory layers]
    dCA1_us_gated <- function(dVTA_val, previousdCA1) {
      if (dVTA_val <= 0) return(0)  # No US surprise → no learning signal for sensory layers

      d <- 0; n <- 0
      for (i in 1:length(network)) {
        if (network[[i]]@Layer == "Hippocampal") {
          d <- d + abs(network[[i]]@Activation - network[[i]]@PreviousActivation)
          n <- n + 1
        }
      }
      dH <- ifelse(n == 0, 0, d / n)
      if (missing(previousdCA1)) previousdCA1 <- 0
      return(dH + dVTA_val * (1 - previousdCA1))
    }

    Compute.r <- function(npe) {
      swe <- 0; swi <- 0
      if (length(network[[npe]]@InputConnections) > 0) {
        for (i in 1:length(network[[npe]]@InputConnections)) {
          pre <- network[[npe]]@InputConnections[[i]]@preSinapticNPE
          if (network[[pre]]@Type == "Excitatory") {
            swe <- swe + ifelse(network[[pre]]@Layer == "US", 0, network[[npe]]@InputConnections[[i]]@weight)
          } else { swi <- swi + network[[npe]]@InputConnections[[i]]@weight }
        }
      }
      return(c(1 - swe, 1 - swi))
    }

    network <- list()
    for (i in 1:nrow(NPEs)) {
      network[[NPEs[i, 1]]] <- new("NPE",
        Name = NPEs[i, 1], Type = NPEs[i, 2], Layer = NPEs[i, 3],
        Activation = NPEs[i, 4], TemporalSummation = NPEs[i, 5],
        ActivationDecay = NPEs[i, 6], mu = NPEs[i, 7],
        sigma = NPEs[i, 8], logisSigma = NPEs[i, 9],
        PreviousExcitatoryInput = 0, ExcitatoryInput = 0,
        InhibitoryInput = 0, PreviousActivation = 0)
    }
    for (i in 1:nrow(Connections)) {
      cn <- paste(Connections[i, 1], Connections[i, 2], sep = "-")
      post <- Connections[i, 2]
      network[[post]]@InputConnections[[cn]] <- new("Connection",
        Name = cn, weight = Connections[i, 3],
        alpha = Connections[i, 4], beta = Connections[i, 5],
        alpha_prime = Connections[i, 6], beta_prime = Connections[i, 7],
        preSinapticNPE = Connections[i, 1])
    }

    n.input <- unlist(lapply(network, function(x) x@Layer %in% c("PrimarySensory", "US")))
    if (!any(names(saveData) == "Elements")) {
      saveData$Elements <- c(
        unname(unlist(lapply(network, FUN = function(x) if (!x@Layer %in% c("PrimarySensory", "US")) return(x@Name)))),
        unname(unlist(lapply(network, FUN = function(x) return(lapply(x@InputConnections, FUN = function(y) return(y@Name)))))))
    }
    if (!any(names(saveData) == "TimeSteps")) saveData$TimeSteps <- unique(TimeSteps$TimeStep)
    saveData$Elements <- union(names(n.input)[n.input], saveData$Elements)

    data.sim <- as.data.frame(matrix(nrow = length(which(TimeSteps$TimeStep %in% saveData$TimeSteps)), ncol = length(saveData$Elements) + 5))
    colnames(data.sim) <- c("Phase", "Trial", "TimeStep", saveData$Elements, "dVTA", "dH")

    LearningRuleIsActive <- F; PreviousdCA1 <- 0; current.phase <- 1; t <- 0

    for (ts in 1:nrow(TimeSteps)) {
      if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) {
        t <- t + 1; data.sim[t, 1] <- TimeSteps[ts, 1]; data.sim[t, 2] <- TimeSteps[ts, 2]; data.sim[t, 3] <- TimeSteps[ts, 3]
      }
      if (ts > 1 && TimeSteps[ts - 1, 1] != TimeSteps[ts, 1]) current.phase <- current.phase + 1
      if (TimeSteps[ts, 3] == 1 & !HasITI[current.phase]) {
        for (i in 1:length(network)) { network[[i]]@Activation <- L(0, network[[i]]@logisSigma); network[[i]]@ExcitatoryInput <- L(0, network[[i]]@logisSigma) }
      }
      for (npu in 1:length(network)) { if (network[[npu]]@Layer == "US" | network[[npu]]@Layer == "PrimarySensory") network[[npu]]@Activation <- 0 }
      LearningRuleIsActive <- as.logical(TimeSteps[ts, ncol(TimeSteps)])
      for (unit in seq(4, ncol(TimeSteps) - 1, 2)) network[[TimeSteps[ts, unit]]]@Activation <- TimeSteps[ts, unit + 1]

      scrambledNPEs <- sample(1:length(network), length(network), replace = F)
      for (i in scrambledNPEs) {
        network[[i]]@PreviousActivation <- network[[i]]@Activation
        network[[i]]@PreviousExcitatoryInput <- network[[i]]@ExcitatoryInput
        if (network[[i]]@Layer == "US" | network[[i]]@Layer == "PrimarySensory") {
          if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) data.sim[t, network[[i]]@Name] <- network[[i]]@Activation; next
        }
        npe.unc <- FALSE; us.act <- 0
        if (network[[i]]@Layer %in% c("Dopaminergic", "PrimaryMotor") && length(network[[i]]@InputConnections) > 0) {
          for (j in 1:length(network[[i]]@InputConnections)) {
            pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
            if (network[[pre]]@Layer == "US" && network[[pre]]@Activation > 0) { npe.unc <- TRUE; us.act <- network[[pre]]@Activation; break }
          }
        }
        if (npe.unc) { network[[i]]@Activation <- us.act } else {
          if (threshold == "gaussian") { network[[i]]@Threshold <- rnorm(1, network[[i]]@mu, network[[i]]@sigma)
          } else { p <- estBetaParams(network[[i]]@mu, network[[i]]@sigma); network[[i]]@Threshold <- rbeta(1, p$alpha, p$beta) }
          inputs <- ComputeInputs(network[[i]]); network[[i]]@ExcitatoryInput <- inputs[1]; network[[i]]@InhibitoryInput <- inputs[2]
          pe <- L(network[[i]]@ExcitatoryInput, network[[i]]@logisSigma); pi <- L(network[[i]]@InhibitoryInput, network[[i]]@logisSigma)
          if (network[[i]]@Layer != "PrimarySensory" & network[[i]]@Layer != "US") {
            if (pe > pi) { if (pe >= network[[i]]@Threshold) {
              network[[i]]@Activation <- pe + network[[i]]@TemporalSummation * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) * (1 - pe) - pi
            } else { network[[i]]@Activation <- L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) - network[[i]]@ActivationDecay * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) }
            } else { network[[i]]@Activation <- 0 }
          }
        }
        if (network[[i]]@Name %in% saveData$Elements & TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) data.sim[t, network[[i]]@Name] <- network[[i]]@Activation
      }
      if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) { data.sim[t, "dVTA"] <- 0; data.sim[t, "dH"] <- 0 }

      if (LearningRuleIsActive) {
        dD <- dVTA_us_gated()
        dH <- dCA1_us_gated(dD, PreviousdCA1); PreviousdCA1 <- dH
        if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) { data.sim[t, "dVTA"] <- dD; data.sim[t, "dH"] <- dH }
        if (pupdate == "async_random" || pupdate == "sync_random") { scrambledNPEs <- sample(1:length(network), length(network), replace = F) } else { scrambledNPEs <- 1:length(network) }

        for (i in scrambledNPEs) {
          inputs <- ComputeInputs(network[[i]]); network[[i]]@ExcitatoryInput <- inputs[1]; network[[i]]@InhibitoryInput <- inputs[2]
          if (network[[i]]@Layer == "US" | network[[i]]@Layer == "PrimarySensory") next
          network[[i]]@r <- Compute.r(i)
          if (pupdate == "async_random" || pupdate == "sync_random") {
            sc <- sample(1:length(network[[i]]@InputConnections), length(network[[i]]@InputConnections), replace = F)
          } else { sc <- 1:length(network[[i]]@InputConnections) }
          if (network[[i]]@Layer %in% c("AssociativeSensory", "Hippocampal")) { d <- dH } else { d <- dD }

          for (j in sc) {
            pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
            if (network[[pre]]@Layer != "US") {
              if (d >= disc) {
                if (network[[pre]]@Type == "Excitatory") {
                  network[[i]]@InputConnections[[j]]@p <- ifelse(network[[i]]@ExcitatoryInput == 0, 0, network[[pre]]@Activation * network[[i]]@InputConnections[[j]]@weight / network[[i]]@ExcitatoryInput)
                  network[[i]]@InputConnections[[j]]@weight <- network[[i]]@InputConnections[[j]]@weight + network[[i]]@InputConnections[[j]]@alpha * network[[i]]@r[1] * network[[i]]@Activation * d * network[[i]]@InputConnections[[j]]@p
                } else {
                  network[[i]]@InputConnections[[j]]@p <- ifelse(network[[i]]@InhibitoryInput == 0, 0, network[[pre]]@Activation * network[[i]]@InputConnections[[j]]@weight / network[[i]]@InhibitoryInput)
                  network[[i]]@InputConnections[[j]]@weight <- network[[i]]@InputConnections[[j]]@weight + network[[i]]@InputConnections[[j]]@alpha_prime * network[[i]]@r[2] * network[[i]]@Activation * d * network[[i]]@InputConnections[[j]]@p
                }
              } else {
                network[[i]]@InputConnections[[j]]@weight <- network[[i]]@InputConnections[[j]]@weight -
                  ifelse(network[[pre]]@Type == "Excitatory", network[[i]]@InputConnections[[j]]@beta, network[[i]]@InputConnections[[j]]@beta_prime) *
                  network[[i]]@InputConnections[[j]]@weight * network[[pre]]@Activation * network[[i]]@Activation
              }
            }
            network[[i]]@InputConnections[[j]]@weight <- min(max(network[[i]]@InputConnections[[j]]@weight, 0), 1)
            if (network[[i]]@InputConnections[[j]]@Name %in% saveData$Elements & TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps)
              data.sim[t, network[[i]]@InputConnections[[j]]@Name] <- network[[i]]@InputConnections[[j]]@weight
          }
        }
      } else {
        for (i in scrambledNPEs) {
          if (pupdate == "async_random" || pupdate == "sync_random") { scs <- sample(1:length(network[[i]]@InputConnections), length(network[[i]]@InputConnections), replace = F) } else { scs <- 1:length(network[[i]]@InputConnections) }
          for (j in scs) { if (network[[i]]@InputConnections[[j]]@Name %in% saveData$Elements & TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) data.sim[t, network[[i]]@InputConnections[[j]]@Name] <- network[[i]]@InputConnections[[j]]@weight }
        }
      }
    }
    return(data.sim)
  }, error = function(e) stop(e))
}

# ==================== TEST ====================
template <- get_template_data("extinction")
n_networks <- 10

cat("=== EXPERIMENT M: Full DTD Adjusted Fix ===\n")
cat("Fix 1: US-gated dVTA (only US surplus counts)\n")
cat("Fix 2: US-gated dCA1 (dCA1=0 when no US)\n\n")

results <- data.frame(M1D_t = numeric(), M1D_e = numeric(), Pct = numeric())

for (net in 1:n_networks) {
  TS <- Create.Phases(phases = template$contingencies, trials = template$trials)
  res <- Simulate.DBP.FullFix(
    NPEs = template$npes, Connections = template$connections,
    TimeSteps = TS, HasITI = template$hasITI, disc = 0.001)
  tr <- which(res$Phase == "training"); ex <- which(res$Phase == "extinction")
  mt <- res[max(tr), "M..1-D"]; me <- res[max(ex), "M..1-D"]
  pct <- (me - mt) / mt * 100
  results <- rbind(results, data.frame(M1D_t = mt, M1D_e = me, Pct = pct))
  cat(sprintf("  Net %d: train=%.4f, ext=%.4f, change=%.1f%%\n", net, mt, me, pct))
}

cat(sprintf("\nMean training:  %.4f (SD=%.4f)\n", mean(results$M1D_t), sd(results$M1D_t)))
cat(sprintf("Mean extinction: %.4f (SD=%.4f)\n", mean(results$M1D_e), sd(results$M1D_e)))
cat(sprintf("Mean change:     %.1f%% (SD=%.1f%%)\n", mean(results$Pct), sd(results$Pct)))
cat(sprintf("Nets >50%% ext:  %d/%d\n", sum(results$Pct < -50), n_networks))
cat(sprintf("Nets >20%% ext:  %d/%d\n", sum(results$Pct < -20), n_networks))

# Show trajectory and all weights for net 1
cat("\n=== Per-trial trajectory (net 1) ===\n")
set.seed(42)
TS <- Create.Phases(phases = template$contingencies, trials = template$trials)
res <- Simulate.DBP.FullFix(NPEs = template$npes, Connections = template$connections,
  TimeSteps = TS, HasITI = template$hasITI, disc = 0.001)
tr <- which(res$Phase == "training"); ex <- which(res$Phase == "extinction")

cat("Phase     Trial  M''1->D\n")
for (trial in c(1, 10, 50, 80, 100)) {
  rows <- tr[res[tr, "Trial"] == trial]
  if (length(rows) > 0) cat(sprintf("Training  %3d    %.4f\n", trial, res[max(rows), "M..1-D"]))
}
cat("---\n")
for (trial in c(1, 5, 10, 20, 30, 50, 70, 100)) {
  rows <- ex[res[ex, "Trial"] == trial]
  if (length(rows) > 0) cat(sprintf("Extinct   %3d    %.4f\n", trial, res[max(rows), "M..1-D"]))
}

cat("\n--- All connection weights ---\n")
cc <- grep("-", colnames(res), value = TRUE)
cat(sprintf("%-20s %12s %12s %12s\n", "Connection", "End_train", "End_ext", "Change"))
for (col in cc) {
  cat(sprintf("%-20s %12.6f %12.6f %12.6f\n", col, res[max(tr), col], res[max(ex), col], res[max(ex), col] - res[max(tr), col]))
}
