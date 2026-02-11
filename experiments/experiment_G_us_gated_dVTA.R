# Experiment G: US-Gated dVTA (DTD Adjusted Fix v2)
#
# Problem: During extinction, D activates via conditioned pathway (M''1->D),
#          creating positive dVTA that triggers INCREMENT.
#
# Fix: Compare D activation to what it would be WITHOUT US input.
#      dVTA = D_actual - D_conditioned
#      where D_conditioned = D computed only from non-US inputs
#      This means dVTA is only positive when the US ACTUALLY contributes
#      something beyond what the CS already predicts.
#
# During training (CS+US present): D_actual=1.0 (US drives D), D_conditioned~0.7
#   => dVTA = 1.0 - 0.7 = +0.3 => INCREMENT (good!)
# During extinction (CS only): D_actual=D_conditioned (no US)
#   => dVTA = 0 => DECREMENT (correct!)
#
# This is exactly the Rescorla-Wagner prediction error:
#   λ (actual US) - ΣV (predicted US via conditioned connections)

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

set.seed(42)

Simulate.DBP.USGated <- function(NPEs, Connections, TimeSteps, HasITI = c(FALSE),
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

    L <- function(x, sigma) {
      return(1 / (1 + exp((-x + 0.5) / sigma)))
    }

    estBetaParams <- function(mu, sigma) {
      var <- sigma^2
      alpha <- mu * (mu - mu^2 - var) / var
      beta <- (mu * (1 - mu) / var - 1) * (1 - mu)
      if (alpha <= 0 | beta <= 0) stop("invalid beta params")
      return(params = list(alpha = alpha, beta = beta))
    }

    # NEW: US-gated dVTA
    # Computes what D WOULD activate as, only counting US-driven input
    dVTA_us_gated <- function() {
      d_actual <- 0
      d_conditioned <- 0
      n <- 0

      for (i in 1:length(network)) {
        if (network[[i]]@Layer == "Dopaminergic") {
          d_actual <- d_actual + network[[i]]@Activation
          n <- n + 1

          # Compute what D activation would be without US
          # Check if this D unit is currently being driven by US
          us_driving <- FALSE
          for (j in 1:length(network[[i]]@InputConnections)) {
            pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
            if (network[[pre]]@Layer == "US" && network[[pre]]@Activation > 0) {
              us_driving <- TRUE
              break
            }
          }

          if (us_driving) {
            # D is being driven by US: compute conditioned-only activation
            # D would have activated to the conditioned pathway level
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
            # Compute what D activation would be
            p_epsp <- L(exc_input_no_us, network[[i]]@logisSigma)
            p_ipsp <- L(inh_input, network[[i]]@logisSigma)
            if (p_epsp > p_ipsp && p_epsp >= network[[i]]@Threshold) {
              d_conditioned <- d_conditioned + p_epsp + network[[i]]@TemporalSummation *
                L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) * (1 - p_epsp) - p_ipsp
            } else if (p_epsp > p_ipsp) {
              d_conditioned <- d_conditioned + L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) -
                network[[i]]@ActivationDecay * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma)
            } else {
              d_conditioned <- d_conditioned + 0
            }
          } else {
            # D is NOT being driven by US, so conditioned = actual
            d_conditioned <- d_conditioned + network[[i]]@Activation
          }
        }
      }
      if (n == 0) return(0)
      # Prediction error: actual - conditioned prediction
      return((d_actual - d_conditioned) / n)
    }

    dCA1 <- function(dVTA, previousdCA1) {
      d <- 0
      n <- 0
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
      sum.weights.exc <- 0
      sum.weights.inh <- 0
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

    # Build network
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
    if (!any(names(saveData) == "TimeSteps")) {
      saveData$TimeSteps <- unique(TimeSteps$TimeStep)
    }
    saveData$Elements <- union(names(n.input)[n.input], saveData$Elements)

    data.sim <- as.data.frame(matrix(nrow = length(which(TimeSteps$TimeStep %in% saveData$TimeSteps)), ncol = length(saveData$Elements) + 5))
    colnames(data.sim) <- c("Phase", "Trial", "TimeStep", saveData$Elements, "dVTA", "dH")

    LearningRuleIsActive <- F
    PreviousdCA1 <- 0
    current.phase <- 1
    t <- 0

    for (ts in 1:nrow(TimeSteps)) {
      if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) {
        t <- t + 1
        data.sim[t, 1] <- TimeSteps[ts, 1]
        data.sim[t, 2] <- TimeSteps[ts, 2]
        data.sim[t, 3] <- TimeSteps[ts, 3]
      }

      if (ts > 1 && TimeSteps[ts - 1, 1] != TimeSteps[ts, 1]) {
        current.phase <- current.phase + 1
      }

      if (TimeSteps[ts, 3] == 1 & !HasITI[current.phase]) {
        for (i in 1:length(network)) {
          network[[i]]@Activation <- L(0, network[[i]]@logisSigma)
          network[[i]]@ExcitatoryInput <- L(0, network[[i]]@logisSigma)
        }
      }

      for (npu in 1:length(network)) {
        if (network[[npu]]@Layer == "US" | network[[npu]]@Layer == "PrimarySensory") {
          network[[npu]]@Activation <- 0
        }
      }

      LearningRuleIsActive <- as.logical(TimeSteps[ts, ncol(TimeSteps)])

      for (unit in seq(4, ncol(TimeSteps) - 1, 2)) {
        network[[TimeSteps[ts, unit]]]@Activation <- TimeSteps[ts, unit + 1]
      }

      scrambledNPEs <- sample(1:length(network), length(network), replace = F)

      for (i in scrambledNPEs) {
        currentPostSinapticNPE <- network[[i]]@Name
        network[[i]]@PreviousActivation <- network[[i]]@Activation
        network[[i]]@PreviousExcitatoryInput <- network[[i]]@ExcitatoryInput

        if (network[[i]]@Layer == "US" | network[[i]]@Layer == "PrimarySensory") {
          if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) data.sim[t, network[[i]]@Name] <- network[[i]]@Activation
          next
        }

        npe.is.unconditionally.activated <- FALSE
        us.activation <- 0

        if (network[[i]]@Layer %in% c("Dopaminergic", "PrimaryMotor") && length(network[[i]]@InputConnections) > 0) {
          for (j in 1:length(network[[i]]@InputConnections)) {
            pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
            if (network[[pre]]@Layer == "US" && network[[pre]]@Activation > 0) {
              npe.is.unconditionally.activated <- TRUE
              us.activation <- network[[pre]]@Activation
              break
            }
          }
        }

        if (npe.is.unconditionally.activated) {
          network[[i]]@Activation <- us.activation
        } else {
          if (threshold == "gaussian") {
            network[[i]]@Threshold <- rnorm(1, network[[i]]@mu, network[[i]]@sigma)
          } else {
            p <- estBetaParams(network[[i]]@mu, network[[i]]@sigma)
            network[[i]]@Threshold <- rbeta(1, p$alpha, p$beta)
          }

          inputs <- ComputeInputs(network[[i]])
          network[[i]]@ExcitatoryInput <- inputs[1]
          network[[i]]@InhibitoryInput <- inputs[2]

          p_epsp <- L(network[[i]]@ExcitatoryInput, network[[i]]@logisSigma)
          p_ipsp <- L(network[[i]]@InhibitoryInput, network[[i]]@logisSigma)

          if (network[[i]]@Layer != "PrimarySensory" & network[[i]]@Layer != "US") {
            if (p_epsp > p_ipsp) {
              if (p_epsp >= network[[i]]@Threshold) {
                network[[i]]@Activation <- p_epsp + network[[i]]@TemporalSummation * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) * (1 - p_epsp) - p_ipsp
              } else {
                network[[i]]@Activation <- L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) - network[[i]]@ActivationDecay * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma)
              }
            } else {
              network[[i]]@Activation <- 0
            }
          }
        }

        if (network[[i]]@Name %in% saveData$Elements & TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) {
          data.sim[t, network[[i]]@Name] <- network[[i]]@Activation
        }
      }

      if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) {
        data.sim[t, "dVTA"] <- 0
        data.sim[t, "dH"] <- 0
      }

      if (LearningRuleIsActive) {
        # Use US-gated dVTA
        dD <- dVTA_us_gated()
        dH <- dCA1(dD, PreviousdCA1)
        PreviousdCA1 <- dH

        if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) {
          data.sim[t, "dVTA"] <- dD
          data.sim[t, "dH"] <- dH
        }

        if (pupdate == "async_random" || pupdate == "sync_random") {
          scrambledNPEs <- sample(1:length(network), length(network), replace = F)
        } else {
          scrambledNPEs <- 1:length(network)
        }

        for (i in scrambledNPEs) {
          currentPostSinapticNPE <- network[[i]]@Name
          inputs <- ComputeInputs(network[[i]])
          network[[i]]@ExcitatoryInput <- inputs[1]
          network[[i]]@InhibitoryInput <- inputs[2]

          if (network[[currentPostSinapticNPE]]@Layer == "US" | network[[currentPostSinapticNPE]]@Layer == "PrimarySensory") {
            next
          }

          network[[currentPostSinapticNPE]]@r <- Compute.r(i)

          if (pupdate == "async_random" || pupdate == "sync_random") {
            scrambledConnections <- sample(1:length(network[[currentPostSinapticNPE]]@InputConnections), length(network[[currentPostSinapticNPE]]@InputConnections), replace = F)
          } else {
            scrambledConnections <- 1:length(network[[currentPostSinapticNPE]]@InputConnections)
          }

          if (network[[currentPostSinapticNPE]]@Layer %in% c("AssociativeSensory", "Hippocampal")) {
            d <- dH
          } else {
            d <- dD
          }

          for (j in scrambledConnections) {
            pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE

            if (network[[pre]]@Layer != "US") {
              if (d >= disc) {
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
                # Keep original exponential decrement
                network[[i]]@InputConnections[[j]]@weight <- network[[i]]@InputConnections[[j]]@weight -
                  ifelse(network[[pre]]@Type == "Excitatory",
                    network[[i]]@InputConnections[[j]]@beta,
                    network[[i]]@InputConnections[[j]]@beta_prime
                  ) *
                    network[[i]]@InputConnections[[j]]@weight *
                    network[[pre]]@Activation *
                    network[[i]]@Activation
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
          } else {
            scrambleConnections <- 1:length(network[[i]]@InputConnections)
          }
          for (j in scrambleConnections) {
            if (network[[i]]@InputConnections[[j]]@Name %in% saveData$Elements &
              TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) {
              data.sim[t, network[[i]]@InputConnections[[j]]@Name] <- network[[i]]@InputConnections[[j]]@weight
            }
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

cat("=== EXPERIMENT G: US-Gated dVTA ===\n")
cat("dVTA = D_actual - D_conditioned (only US-driven surplus counts)\n\n")

results <- data.frame(
  Network = integer(),
  M1D_end_training = numeric(),
  M1D_end_extinction = numeric(),
  M1D_change = numeric(),
  Pct_change = numeric(),
  stringsAsFactors = FALSE
)

for (net in 1:n_networks) {
  TimeSteps <- Create.Phases(
    phases = template$contingencies,
    trials = template$trials
  )

  result <- Simulate.DBP.USGated(
    NPEs = template$npes,
    Connections = template$connections,
    TimeSteps = TimeSteps,
    HasITI = template$hasITI,
    disc = 0.001
  )

  training_rows <- which(result$Phase == "training")
  extinction_rows <- which(result$Phase == "extinction")

  m1d_col <- "M..1-D"
  m1d_end_training <- result[max(training_rows), m1d_col]
  m1d_end_extinction <- result[max(extinction_rows), m1d_col]
  m1d_change <- m1d_end_extinction - m1d_end_training
  pct_change <- m1d_change / m1d_end_training * 100

  results <- rbind(results, data.frame(
    Network = net,
    M1D_end_training = round(m1d_end_training, 6),
    M1D_end_extinction = round(m1d_end_extinction, 6),
    M1D_change = round(m1d_change, 6),
    Pct_change = round(pct_change, 2),
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  Net %d: M''1->D train=%.4f, ext=%.4f, change=%.4f (%.1f%%)\n",
              net, m1d_end_training, m1d_end_extinction, m1d_change, pct_change))
}

cat(sprintf("\nMean M''1->D training:  %.4f (SD=%.4f)\n", mean(results$M1D_end_training), sd(results$M1D_end_training)))
cat(sprintf("Mean M''1->D extinction: %.4f (SD=%.4f)\n", mean(results$M1D_end_extinction), sd(results$M1D_end_extinction)))
cat(sprintf("Mean change:            %.4f (%.1f%%)\n", mean(results$M1D_change), mean(results$Pct_change)))
cat(sprintf("Nets >50%% extinction:   %d/%d\n", sum(results$Pct_change < -50), n_networks))
cat(sprintf("Nets >20%% extinction:   %d/%d\n", sum(results$Pct_change < -20), n_networks))

# Show per-trial trajectory for net 1
cat("\n=== Per-trial trajectory (net 1) ===\n")
set.seed(42)
TimeSteps <- Create.Phases(phases = template$contingencies, trials = template$trials)
result <- Simulate.DBP.USGated(
  NPEs = template$npes, Connections = template$connections,
  TimeSteps = TimeSteps, HasITI = template$hasITI, disc = 0.001
)

training_rows <- which(result$Phase == "training")
extinction_rows <- which(result$Phase == "extinction")
m1d_col <- "M..1-D"

cat("Phase     Trial  M''1->D   dVTA\n")
for (trial in c(1, 5, 10, 20, 50, 80, 100)) {
  rows <- training_rows[result[training_rows, "Trial"] == trial]
  if (length(rows) > 0) {
    last_row <- max(rows)
    cat(sprintf("Training  %3d    %.4f   %.4f\n", trial, result[last_row, m1d_col], result[last_row, "dVTA"]))
  }
}
cat("---\n")
for (trial in c(1, 2, 3, 5, 10, 20, 30, 50, 70, 100)) {
  rows <- extinction_rows[result[extinction_rows, "Trial"] == trial]
  if (length(rows) > 0) {
    last_row <- max(rows)
    cat(sprintf("Extinct   %3d    %.4f   %.4f\n", trial, result[last_row, m1d_col], result[last_row, "dVTA"]))
  }
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
