# Experiment H: US-Gated dVTA + Parameter Tuning
#
# Exp G showed that US-gated dVTA correctly:
# - Signals dVTA > 0 only when US is present (training)
# - Signals dVTA = 0 during extinction (no false increments!)
#
# But extinction is too slow (-4.1%) because:
# - DECREMENT fires every timestep during extinction (dVTA=0 < disc)
# - But β=0.12 * w * a_pre * a_post is very small per timestep
#
# Solutions to try:
# 1. Higher β (0.2, 0.5, 1.0)
# 2. Remove w multiplier from decrement (linear decay)
# 3. Both

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

set.seed(42)

Simulate.DBP.Adjusted2 <- function(NPEs, Connections, TimeSteps, HasITI = c(FALSE),
                                    disc = 0.001, threshold = "gaussian",
                                    pupdate = "async_random",
                                    saveData = list(),
                                    beta_multiplier = 1.0,
                                    use_linear_decrement = FALSE) {

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
      var <- sigma^2
      alpha <- mu * (mu - mu^2 - var) / var
      beta <- (mu * (1 - mu) / var - 1) * (1 - mu)
      if (alpha <= 0 | beta <= 0) stop("invalid beta params")
      return(params = list(alpha = alpha, beta = beta))
    }

    # US-gated dVTA
    dVTA_us_gated <- function() {
      d_actual <- 0
      d_conditioned <- 0
      n <- 0
      for (i in 1:length(network)) {
        if (network[[i]]@Layer == "Dopaminergic") {
          d_actual <- d_actual + network[[i]]@Activation
          n <- n + 1
          us_driving <- FALSE
          for (j in 1:length(network[[i]]@InputConnections)) {
            pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
            if (network[[pre]]@Layer == "US" && network[[pre]]@Activation > 0) {
              us_driving <- TRUE
              break
            }
          }
          if (us_driving) {
            exc_input_no_us <- 0
            inh_input <- 0
            for (j in 1:length(network[[i]]@InputConnections)) {
              pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
              if (network[[pre]]@Layer != "US") {
                if (network[[pre]]@Type == "Excitatory") {
                  exc_input_no_us <- exc_input_no_us + network[[pre]]@Activation * network[[i]]@InputConnections[[j]]@weight
                } else {
                  inh_input <- inh_input + network[[pre]]@Activation * network[[i]]@InputConnections[[j]]@weight
                }
              }
            }
            p_epsp <- L(exc_input_no_us, network[[i]]@logisSigma)
            p_ipsp <- L(inh_input, network[[i]]@logisSigma)
            if (p_epsp > p_ipsp && p_epsp >= network[[i]]@Threshold) {
              d_conditioned <- d_conditioned + p_epsp + network[[i]]@TemporalSummation *
                L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) * (1 - p_epsp) - p_ipsp
            } else if (p_epsp > p_ipsp) {
              d_conditioned <- d_conditioned + L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) -
                network[[i]]@ActivationDecay * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma)
            }
          } else {
            d_conditioned <- d_conditioned + network[[i]]@Activation
          }
        }
      }
      if (n == 0) return(0)
      return((d_actual - d_conditioned) / n)
    }

    dCA1 <- function(dVTA, previousdCA1) {
      d <- 0; n <- 0
      for (i in 1:length(network)) {
        if (network[[i]]@Layer == "Hippocampal") {
          d <- d + abs(network[[i]]@Activation - network[[i]]@PreviousActivation)
          n <- n + 1
        }
      }
      dH <- ifelse(n == 0, 0, d / n)
      if (missing(previousdCA1)) previousdCA1 <- 0
      return(dH + dVTA * (1 - previousdCA1))
    }

    Compute.r <- function(npe) {
      sum.weights.exc <- 0; sum.weights.inh <- 0
      if (length(network[[npe]]@InputConnections) > 0) {
        for (i in 1:length(network[[npe]]@InputConnections)) {
          pre <- network[[npe]]@InputConnections[[i]]@preSinapticNPE
          if (network[[pre]]@Type == "Excitatory") {
            sum.weights.exc <- sum.weights.exc + ifelse(network[[pre]]@Layer == "US", 0, network[[npe]]@InputConnections[[i]]@weight)
          } else {
            sum.weights.inh <- sum.weights.inh + network[[npe]]@InputConnections[[i]]@weight
          }
        }
      }
      return(c(1 - sum.weights.exc, 1 - sum.weights.inh))
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
      connection.name <- paste(Connections[i, 1], Connections[i, 2], sep = "-")
      currentPostSinapticNPE <- Connections[i, 2]
      network[[currentPostSinapticNPE]]@InputConnections[[connection.name]] <- new("Connection",
        Name = connection.name, weight = Connections[i, 3],
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
        for (i in 1:length(network)) {
          network[[i]]@Activation <- L(0, network[[i]]@logisSigma)
          network[[i]]@ExcitatoryInput <- L(0, network[[i]]@logisSigma)
        }
      }
      for (npu in 1:length(network)) {
        if (network[[npu]]@Layer == "US" | network[[npu]]@Layer == "PrimarySensory") network[[npu]]@Activation <- 0
      }
      LearningRuleIsActive <- as.logical(TimeSteps[ts, ncol(TimeSteps)])
      for (unit in seq(4, ncol(TimeSteps) - 1, 2)) {
        network[[TimeSteps[ts, unit]]]@Activation <- TimeSteps[ts, unit + 1]
      }
      scrambledNPEs <- sample(1:length(network), length(network), replace = F)
      for (i in scrambledNPEs) {
        network[[i]]@PreviousActivation <- network[[i]]@Activation
        network[[i]]@PreviousExcitatoryInput <- network[[i]]@ExcitatoryInput
        if (network[[i]]@Layer == "US" | network[[i]]@Layer == "PrimarySensory") {
          if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) data.sim[t, network[[i]]@Name] <- network[[i]]@Activation
          next
        }
        npe.is.unconditionally.activated <- FALSE; us.activation <- 0
        if (network[[i]]@Layer %in% c("Dopaminergic", "PrimaryMotor") && length(network[[i]]@InputConnections) > 0) {
          for (j in 1:length(network[[i]]@InputConnections)) {
            pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
            if (network[[pre]]@Layer == "US" && network[[pre]]@Activation > 0) {
              npe.is.unconditionally.activated <- TRUE; us.activation <- network[[pre]]@Activation; break
            }
          }
        }
        if (npe.is.unconditionally.activated) {
          network[[i]]@Activation <- us.activation
        } else {
          if (threshold == "gaussian") { network[[i]]@Threshold <- rnorm(1, network[[i]]@mu, network[[i]]@sigma)
          } else { p <- estBetaParams(network[[i]]@mu, network[[i]]@sigma); network[[i]]@Threshold <- rbeta(1, p$alpha, p$beta) }
          inputs <- ComputeInputs(network[[i]])
          network[[i]]@ExcitatoryInput <- inputs[1]; network[[i]]@InhibitoryInput <- inputs[2]
          p_epsp <- L(network[[i]]@ExcitatoryInput, network[[i]]@logisSigma)
          p_ipsp <- L(network[[i]]@InhibitoryInput, network[[i]]@logisSigma)
          if (network[[i]]@Layer != "PrimarySensory" & network[[i]]@Layer != "US") {
            if (p_epsp > p_ipsp) {
              if (p_epsp >= network[[i]]@Threshold) {
                network[[i]]@Activation <- p_epsp + network[[i]]@TemporalSummation * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) * (1 - p_epsp) - p_ipsp
              } else {
                network[[i]]@Activation <- L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) - network[[i]]@ActivationDecay * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma)
              }
            } else { network[[i]]@Activation <- 0 }
          }
        }
        if (network[[i]]@Name %in% saveData$Elements & TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps)
          data.sim[t, network[[i]]@Name] <- network[[i]]@Activation
      }
      if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) { data.sim[t, "dVTA"] <- 0; data.sim[t, "dH"] <- 0 }

      if (LearningRuleIsActive) {
        dD <- dVTA_us_gated()
        dH <- dCA1(dD, PreviousdCA1); PreviousdCA1 <- dH
        if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) { data.sim[t, "dVTA"] <- dD; data.sim[t, "dH"] <- dH }

        if (pupdate == "async_random" || pupdate == "sync_random") {
          scrambledNPEs <- sample(1:length(network), length(network), replace = F)
        } else { scrambledNPEs <- 1:length(network) }

        for (i in scrambledNPEs) {
          inputs <- ComputeInputs(network[[i]])
          network[[i]]@ExcitatoryInput <- inputs[1]; network[[i]]@InhibitoryInput <- inputs[2]
          if (network[[i]]@Layer == "US" | network[[i]]@Layer == "PrimarySensory") next
          network[[i]]@r <- Compute.r(i)
          if (pupdate == "async_random" || pupdate == "sync_random") {
            scrambledConnections <- sample(1:length(network[[i]]@InputConnections), length(network[[i]]@InputConnections), replace = F)
          } else { scrambledConnections <- 1:length(network[[i]]@InputConnections) }
          if (network[[i]]@Layer %in% c("AssociativeSensory", "Hippocampal")) { d <- dH } else { d <- dD }

          for (j in scrambledConnections) {
            pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
            if (network[[pre]]@Layer != "US") {
              if (d >= disc) {
                # INCREMENT (same as original)
                if (network[[pre]]@Type == "Excitatory") {
                  network[[i]]@InputConnections[[j]]@p <- ifelse(network[[i]]@ExcitatoryInput == 0, 0, network[[pre]]@Activation * network[[i]]@InputConnections[[j]]@weight / network[[i]]@ExcitatoryInput)
                  network[[i]]@InputConnections[[j]]@weight <- network[[i]]@InputConnections[[j]]@weight +
                    network[[i]]@InputConnections[[j]]@alpha * network[[i]]@r[1] *
                      network[[i]]@Activation * d * network[[i]]@InputConnections[[j]]@p
                } else {
                  network[[i]]@InputConnections[[j]]@p <- ifelse(network[[i]]@InhibitoryInput == 0, 0, network[[pre]]@Activation * network[[i]]@InputConnections[[j]]@weight / network[[i]]@InhibitoryInput)
                  network[[i]]@InputConnections[[j]]@weight <- network[[i]]@InputConnections[[j]]@weight +
                    network[[i]]@InputConnections[[j]]@alpha_prime * network[[i]]@r[2] *
                      network[[i]]@Activation * d * network[[i]]@InputConnections[[j]]@p
                }
              } else {
                # DECREMENT - may be exponential or linear, with beta_multiplier
                beta_val <- ifelse(network[[pre]]@Type == "Excitatory",
                  network[[i]]@InputConnections[[j]]@beta,
                  network[[i]]@InputConnections[[j]]@beta_prime
                ) * beta_multiplier

                if (use_linear_decrement) {
                  # Linear: Δw = -β * a_pre * a_post (no w multiplier)
                  network[[i]]@InputConnections[[j]]@weight <- network[[i]]@InputConnections[[j]]@weight -
                    beta_val * network[[pre]]@Activation * network[[i]]@Activation
                } else {
                  # Exponential: Δw = -β * w * a_pre * a_post (original)
                  network[[i]]@InputConnections[[j]]@weight <- network[[i]]@InputConnections[[j]]@weight -
                    beta_val * network[[i]]@InputConnections[[j]]@weight *
                    network[[pre]]@Activation * network[[i]]@Activation
                }
              }
            }
            network[[i]]@InputConnections[[j]]@weight <- min(max(network[[i]]@InputConnections[[j]]@weight, 0), 1)
            if (network[[i]]@InputConnections[[j]]@Name %in% saveData$Elements &
              TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) {
              data.sim[t, network[[i]]@InputConnections[[j]]@Name] <- network[[i]]@InputConnections[[j]]@weight
            }
          }
        }
      } else {
        for (i in scrambledNPEs) {
          if (pupdate == "async_random" || pupdate == "sync_random") {
            scrambleConnections <- sample(1:length(network[[i]]@InputConnections), length(network[[i]]@InputConnections), replace = F)
          } else { scrambleConnections <- 1:length(network[[i]]@InputConnections) }
          for (j in scrambleConnections) {
            if (network[[i]]@InputConnections[[j]]@Name %in% saveData$Elements &
              TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps)
              data.sim[t, network[[i]]@InputConnections[[j]]@Name] <- network[[i]]@InputConnections[[j]]@weight
          }
        }
      }
    }
    return(data.sim)
  }, error = function(e) stop(e))
}

# ==================== TEST ====================
template <- get_template_data("extinction")
n_networks <- 10

cat("=== EXPERIMENT H: US-Gated dVTA + Parameter Tuning ===\n\n")

conditions <- list(
  list(name = "US-gated + exp β=0.12 (original)", beta_mult = 1.0, linear = FALSE),
  list(name = "US-gated + exp β=0.5", beta_mult = 0.5/0.12, linear = FALSE),
  list(name = "US-gated + exp β=1.0", beta_mult = 1.0/0.12, linear = FALSE),
  list(name = "US-gated + linear β=0.12", beta_mult = 1.0, linear = TRUE),
  list(name = "US-gated + linear β=0.05", beta_mult = 0.05/0.12, linear = TRUE),
  list(name = "US-gated + linear β=0.02", beta_mult = 0.02/0.12, linear = TRUE)
)

for (cond in conditions) {
  cat(sprintf("--- %s ---\n", cond$name))

  results <- data.frame(
    M1D_end_training = numeric(),
    M1D_end_extinction = numeric(),
    Pct_change = numeric(),
    stringsAsFactors = FALSE
  )

  for (net in 1:n_networks) {
    TimeSteps <- Create.Phases(phases = template$contingencies, trials = template$trials)
    result <- Simulate.DBP.Adjusted2(
      NPEs = template$npes, Connections = template$connections,
      TimeSteps = TimeSteps, HasITI = template$hasITI,
      disc = 0.001, beta_multiplier = cond$beta_mult, use_linear_decrement = cond$linear
    )
    training_rows <- which(result$Phase == "training")
    extinction_rows <- which(result$Phase == "extinction")
    m1d_col <- "M..1-D"
    m1d_t <- result[max(training_rows), m1d_col]
    m1d_e <- result[max(extinction_rows), m1d_col]
    results <- rbind(results, data.frame(
      M1D_end_training = m1d_t, M1D_end_extinction = m1d_e,
      Pct_change = (m1d_e - m1d_t) / m1d_t * 100))
  }

  cat(sprintf("  Training: %.4f (SD=%.4f)  Extinction: %.4f (SD=%.4f)  Change: %.1f%% (SD=%.1f%%)\n",
    mean(results$M1D_end_training), sd(results$M1D_end_training),
    mean(results$M1D_end_extinction), sd(results$M1D_end_extinction),
    mean(results$Pct_change), sd(results$Pct_change)))
  cat(sprintf("  Nets >50%% ext: %d/%d  Nets >20%% ext: %d/%d\n\n",
    sum(results$Pct_change < -50), n_networks,
    sum(results$Pct_change < -20), n_networks))
}
