#################################################################################################
# Librerias necesarias
library(shiny)
library(shinydashboard)
library(DT)
library(igraph)
library(ggraph)
library(readr)
library(tools)
library(tidygraph)
library(dplyr)
library(stringr)
library(visNetwork)
library(jsonlite)
library(shinyWidgets)
library(shinyjs)
library(plotly)
library(tidyr)
library(shinyBS)


#################################################################################################
### 2024/09/11 # added the possibility of a beta distributed threshold
### 2024/09/12 # added the possibility that a D unit be presynaptic, but just like US efferent connections, D efferent connections do not compete
### 2024/09/16 # added the possibility to change the discrepancy criterion in the function call
### 2026/02/07 # fixed simulation inconsistencies and restored historical beta defaults (beta = beta' = 0.12)
Simulate.DBP <- function(NPEs,
                         Connections,
                         TimeSteps,
                         HasITI,
                         threshold = "gaussian",
                         saveData = list(),
                         disc = 0.001) {
  tryCatch(
    {
      setClass("NPE", slots = c(
        Activation = "numeric",
        PreviousActivation = "numeric",
        ActivationDecay = "numeric",
        ExcitatoryInput = "numeric",
        InhibitoryInput = "numeric",
        PreviousExcitatoryInput = "numeric",
        TemporalSummation = "numeric",
        Name = "character",
        Type = "character",
        Layer = "character",
        Threshold = "numeric",
        mu = "numeric",
        sigma = "numeric",
        logisSigma = "numeric",
        InputConnections = "list",
        r = "numeric"
      ))

      setClass("Connection", slots = c(
        weight = "numeric",
        Name = "character",
        alpha = "numeric",
        beta = "numeric",
        alpha_prime = "numeric",
        beta_prime = "numeric",
        preSinapticNPE = "character",
        p = "numeric"
      ))

      # Accept data.frame or tibble inputs from UI/imports.
      NPEs <- as.data.frame(NPEs, stringsAsFactors = FALSE)
      Connections <- as.data.frame(Connections, stringsAsFactors = FALSE)

      if (ncol(NPEs) < 9) {
        stop("NPEs must include at least 9 columns: NPE, Type, Layer, Activation, Temporal.Summation, Activation.Decay, mu, sigma, logisSigma")
      }
      if (ncol(Connections) < 7) {
        stop("Connections must include at least 7 columns: PreSinapticNPE, PostSinapticNPE, Weight, alpha, beta, alpha_prime, beta_prime")
      }

      NPEs <- NPEs[, 1:9, drop = FALSE]
      Connections <- Connections[, 1:7, drop = FALSE]

      NPEs[, 1:3] <- lapply(NPEs[, 1:3, drop = FALSE], as.character)
      NPEs[, 4:9] <- lapply(NPEs[, 4:9, drop = FALSE], function(x) suppressWarnings(as.numeric(x)))

      Connections[, 1:2] <- lapply(Connections[, 1:2, drop = FALSE], as.character)
      Connections[, 3:7] <- lapply(Connections[, 3:7, drop = FALSE], function(x) suppressWarnings(as.numeric(x)))

      if (any(is.na(NPEs[, 4:9]))) {
        stop("NPE parameters contain non-numeric values in Activation/Temporal.Summation/Activation.Decay/mu/sigma/logisSigma")
      }
      if (any(is.na(Connections[, 3:7]))) {
        stop("Connection parameters contain non-numeric values in Weight/alpha/beta/alpha_prime/beta_prime")
      }


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

        if (alpha <= 0 | beta <= 0) stop("this combination of mean and standard deviation for the threshold results in invalid parameters for a beta distribution")

        return(params = list(alpha = alpha, beta = beta))
      }

      dVTA <- function() {
        d <- 0
        n <- 0
        for (i in 1:length(network)) {
          if (network[[i]]@Layer == "Dopaminergic") {
            d <- d + (network[[i]]@Activation - network[[i]]@PreviousActivation)
            n <- n + 1
          }
        }

        if (n == 0) {
          return(0)
        }
        return(d / n)
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

        if (missing(previousdCA1)) {
          previousdCA1 <- 0
        }

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


      # Initialize network ----
      network <- list()

      for (i in 1:nrow(NPEs)) {
        network[[NPEs[i, 1]]] <- new("NPE",
          Name = NPEs[i, 1],
          Type = NPEs[i, 2],
          Layer = NPEs[i, 3],
          Activation = NPEs[i, 4],
          TemporalSummation = NPEs[i, 5],
          ActivationDecay = NPEs[i, 6],
          mu = NPEs[i, 7],
          sigma = NPEs[i, 8],
          logisSigma = NPEs[i, 9],
          PreviousExcitatoryInput = 0,
          ExcitatoryInput = 0,
          InhibitoryInput = 0,
          PreviousActivation = 0
        )
      }

      for (i in 1:nrow(Connections)) {
        connection.name <- paste(Connections[i, 1], Connections[i, 2], sep = "-")
        currentPostSinapticNPE <- Connections[i, 2]
        network[[currentPostSinapticNPE]]@InputConnections[[connection.name]] <- new("Connection",
          Name = connection.name,
          weight = Connections[i, 3],
          alpha = Connections[i, 4],
          beta = 0.12, # historical default for beta
          alpha_prime = Connections[i, 6],
          beta_prime = 0.12, # historical default for beta prime
          preSinapticNPE = Connections[i, 1]
        )
      }

      n.input <- unlist(lapply(network, function(x) x@Layer %in% c("PrimarySensory", "US")))


      if (class(saveData) != "list") stop("saveData must be a list with at least one member: either (Elements), a vector with npes and connections that you want to be saved, or (TimeSteps), a vector of timesteps to be saved")

      if (!any(names(saveData) == "Elements")) {
        saveData$Elements <- c(
          unname(unlist(lapply(network, FUN = function(x) if (!x@Layer %in% c("PrimarySensory", "US")) {
            return(x@Name)
          }))),
          unname(unlist(lapply(network, FUN = function(x) {
            return(lapply(x@InputConnections, FUN = function(y) {
              return(y@Name)
            }))
          })))
        )
      }

      if (!any(names(saveData) == "TimeSteps")) {
        saveData$TimeSteps <- unique(TimeSteps$TimeStep)
      }

      if (any(!saveData$Elements %in% union(NPEs[, 1], paste(Connections[, 1], Connections[, 2], sep = "-")))) stop("SaveData error. Trying to save at least one NPE or one Connection not in the network")
      if (all(!saveData$TimeSteps %in% unique(TimeSteps$TimeStep))) stop("timesteps to be saved not equal to any created timesteps")


      saveData$Elements <- union(names(n.input)[n.input], saveData$Elements)

      data.sim <- as.data.frame(matrix(nrow = length(which(TimeSteps$TimeStep %in% saveData$TimeSteps)), ncol = length(saveData$Elements) + 3))

      colnames(data.sim) <- c("Phase", "Trial", "TimeStep", saveData$Elements)


      LearningRuleIsActive <- F
      PreviousdCA1 <- 0

      pb <- txtProgressBar(max = nrow(TimeSteps), style = 3)

      current.phase <- 1
      t <- 0

      for (ts in 1:nrow(TimeSteps)) {
        ### save phase, trial and timestep
        if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) {
          t <- t + 1
          data.sim[t, 1] <- TimeSteps[ts, 1]
          data.sim[t, 2] <- TimeSteps[ts, 2]
          data.sim[t, 3] <- TimeSteps[ts, 3]
        }

        if (ts > 1 && TimeSteps[ts - 1, 1] != TimeSteps[ts, 1]) {
          current.phase <- current.phase + 1
        }

        # reset activations ----
        # if HasITI[current.phase] is false and trial onset


        if (TimeSteps[ts, 3] == 1 & !HasITI[current.phase]) {
          for (i in 1:length(network)) {
            network[[i]]@Activation <- L(0, network[[i]]@logisSigma)
            network[[i]]@ExcitatoryInput <- L(0, network[[i]]@logisSigma)
          }
        }

        # reset input activations to zero ----

        for (npu in 1:length(network)) {
          if (network[[npu]]@Layer == "US" | network[[npu]]@Layer == "PrimarySensory") {
            network[[npu]]@Activation <- 0
          }
        }

        # set inputs----

        LearningRuleIsActive <- as.logical(TimeSteps[ts, ncol(TimeSteps)])

        for (unit in seq(4, ncol(TimeSteps) - 1, 2)) {
          network[[TimeSteps[ts, unit]]]@Activation <- TimeSteps[ts, unit + 1]
        }


        # update activations ----
        # asynchronous random

        scrambledNPEs <- sample(1:length(network), length(network), replace = F)

        # scrambledNPEs <- 1:length(network) # asynchronous sequential

        for (i in scrambledNPEs) {
          currentPostSinapticNPE <- network[[i]]@Name

          network[[i]]@PreviousActivation <- network[[i]]@Activation

          network[[i]]@PreviousExcitatoryInput <- network[[i]]@ExcitatoryInput

          if (network[[i]]@Layer == "US" | network[[i]]@Layer == "PrimarySensory") {
            if (TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) data.sim[t, network[[i]]@Name] <- network[[i]]@Activation
            next
          }

          # unconditional activation from active US to D or M' units
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
            # unconditional activation
            network[[i]]@Activation <- us.activation
          } else {
            # conditional activation

            if (threshold == "gaussian") {
              network[[i]]@Threshold <- rnorm(1, network[[i]]@mu, network[[i]]@sigma)
            } else {
              p <- estBetaParams(network[[i]]@mu, network[[i]]@sigma)

              network[[i]]@Threshold <- rbeta(1, p$alpha, p$beta)
            }


            # compute excitatory and inhibitory inputs

            inputs <- ComputeInputs(network[[i]])

            network[[i]]@ExcitatoryInput <- inputs[1]
            network[[i]]@InhibitoryInput <- inputs[2]

            p_epsp <- L(network[[i]]@ExcitatoryInput, network[[i]]@logisSigma)

            p_ipsp <- L(network[[i]]@InhibitoryInput, network[[i]]@logisSigma)

            if (network[[i]]@Layer != "PrimarySensory" & network[[i]]@Layer != "US") {
              if (p_epsp > p_ipsp) {
                if (p_epsp >= network[[i]]@Threshold) {
                  # Reactivation
                  network[[i]]@Activation <- p_epsp + network[[i]]@TemporalSummation * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) * (1 - p_epsp) - p_ipsp
                  # network[[i]]@Activation = L(network[[i]]@Activation + network[[i]]@ExcitatoryInput - network[[i]]@InhibitoryInput,network[[i]]@logisSigma)
                } else {
                  # decay

                  network[[i]]@Activation <- L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) - network[[i]]@ActivationDecay * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma)
                  # network[[i]]@Activation = network[[i]]@Activation - network[[i]]@ActivationDecay * network[[i]]@Activation * (1-network[[i]]@Activation)
                  # network[[i]]@Activation = (1-network[[i]]@ActivationDecay) * network[[i]]@Activation
                }
              } else {
                # inhibition

                network[[i]]@Activation <- 0
              }
            }
          }

          if (network[[i]]@Name %in% saveData$Elements &
            TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) {
            data.sim[t, network[[i]]@Name] <- network[[i]]@Activation
          }
        }

        # update weights ----

        if (LearningRuleIsActive) {
          dD <- dVTA()

          dH <- dCA1(dD, PreviousdCA1)

          PreviousdCA1 <- dH

          scrambledNPEs <- sample(1:length(network), length(network), replace = F)

          for (i in scrambledNPEs) {
            currentPostSinapticNPE <- network[[i]]@Name


            # compute excitatory and inhibitory inputs
            inputs <- ComputeInputs(network[[i]])

            network[[i]]@ExcitatoryInput <- inputs[1]
            network[[i]]@InhibitoryInput <- inputs[2]

            if (network[[currentPostSinapticNPE]]@Layer == "US" | network[[currentPostSinapticNPE]]@Layer == "PrimarySensory") {
              next
            }

            network[[currentPostSinapticNPE]]@r <- Compute.r(i)

            scrambledConnections <- sample(1:length(network[[currentPostSinapticNPE]]@InputConnections), length(network[[currentPostSinapticNPE]]@InputConnections), replace = F)

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
        } else { # learning rule is not active

          for (i in scrambledNPEs) {
            scrambleConnections <- sample(1:length(network[[i]]@InputConnections), length(network[[i]]@InputConnections), replace = F)

            for (j in scrambleConnections) {
              if (network[[i]]@InputConnections[[j]]@Name %in% saveData$Elements &
                TimeSteps[ts, "TimeStep"] %in% saveData$TimeSteps) {
                data.sim[t, network[[i]]@InputConnections[[j]]@Name] <- network[[i]]@InputConnections[[j]]@weight
              }
            }
          }
        }

        setTxtProgressBar(pb, ts)
      }
      return(data.sim)
    },
    error = function(e) stop(e)
  )
}


Create.Phases <- function(phases, trials) {
  # phases is a character vector with comma delimited characters. For example,
  # "Training 1, Random,A+/X-,100-100,False,30,30,ITI entrenamiento"
  # "Training 2, In bulk,AX+,100,False,30,30,ITI entrenamiento"
  # "Test,       In bulk,Prueba X,25,False,30,30,ITI entrenamiento"

  # trials is a list with named elements
  # each element is a named vector of comma delimited characters
  # for example, trials[["A+"]] might be:

  #   [1] "S.1,1,S.2,0.65,S.3,0,US,0,True"
  #   [2] "S.1,1,S.2,0.65,S.3,0,US,0,True"
  #   [3] "S.1,1,S.2,0.65,S.3,0,US,0,True"
  #   [4] "S.1,1,S.2,0.65,S.3,0,US,0,True"
  #   [5] "S.1,1,S.2,0.65,S.3,0,US,1,True"

  timesteps <- as.data.frame(matrix(nrow = 0, ncol = 3 + length(unlist(strsplit(trials[[1]][1], ",")))))

  for (i in 1:length(phases)) {
    current.ts <- 1
    current.trial <- 1

    current.phase <- trimws(unlist(strsplit(phases[i], ",")))

    phase.name <- current.phase[1]
    trial.order <- tolower(current.phase[2])
    trial.types <- trimws(unlist(strsplit(current.phase[3], "/")))

    if (any(!(trial.types %in% names(trials)))) stop("One of more trial names do not match trial names in phases")

    trial.numbers <- as.integer(trimws(unlist(strsplit(current.phase[4], "-"))))
    has.iti <- as.logical(trimws(current.phase[5]))

    if (has.iti) {
      min.ITI <- as.integer(trimws(current.phase[6]))
      max.ITI <- as.integer(trimws(current.phase[7]))
      ITITimestep <- trimws(current.phase[8])
    }

    if (trial.order == "in bulk") {
      for (j in 1:length(trial.types)) {
        for (k in 1:trial.numbers[j]) {
          if (has.iti) {
            if (!(ITITimestep %in% names(trials))) {
              stop("ITI time steps not in trials")
            }

            current.ITI <- ifelse(min.ITI == max.ITI, min.ITI, sample(min.ITI:max.ITI, 1))

            for (ts in 1:current.ITI) {
              timesteps <- rbind(timesteps, c(phase.name, current.trial, current.ts, trimws(unlist(strsplit(trials[[ITITimestep]], ",")))))
              current.ts <- current.ts + 1
            }
          } else {
            current.ts <- 1
          }

          for (ts in 1:length(trials[[trial.types[j]]])) {
            timesteps <- rbind(timesteps, c(phase.name, current.trial, current.ts, trimws(unlist(strsplit(trials[[trial.types[j]]][ts], ",")))))
            current.ts <- current.ts + 1
          }
          current.trial <- current.trial + 1
          current.ts <- 1
        }
      }
    } else if (trial.order == "alternated") {
      for (j in 1:trial.numbers[1]) {
        for (k in 1:length(trial.types)) {
          if (has.iti) {
            if (!(ITITimestep %in% names(trials))) {
              stop("ITI time steps not in trials")
            }

            current.ITI <- ifelse(min.ITI == max.ITI, min.ITI, sample(min.ITI:max.ITI, 1))

            for (ts in 1:current.ITI) {
              timesteps <- rbind(timesteps, c(phase.name, current.trial, current.ts, trimws(unlist(strsplit(trials[[ITITimestep]], ",")))))
              current.ts <- current.ts + 1
            }
          } else {
            current.ts <- 1
          }

          for (ts in 1:length(trials[[trial.types[k]]])) {
            timesteps <- rbind(timesteps, c(phase.name, current.trial, current.ts, trimws(unlist(strsplit(trials[[trial.types[k]]][ts], ",")))))
            current.ts <- current.ts + 1
          }
          current.trial <- current.trial + 1
          current.ts <- 1
        }
      }
    } else if (trial.order == "random") {
      trial.sequence <- vector(mode = "integer")

      for (n in 1:length(trial.numbers)) {
        trial.sequence <- c(trial.sequence, rep(n, trial.numbers[n]))
      }

      trial.sequence <- sample(trial.sequence, length(trial.sequence), replace = F)

      for (j in trial.sequence) {
        if (has.iti) {
          if (!(ITITimestep %in% names(trials))) {
            stop("ITI time steps not in trials")
          }

          current.ITI <- ifelse(min.ITI == max.ITI, min.ITI, sample(min.ITI:max.ITI, 1))

          for (ts in 1:current.ITI) {
            timesteps <- rbind(timesteps, c(phase.name, current.trial, current.ts, trimws(unlist(strsplit(trials[[ITITimestep]], ",")))))
            current.ts <- current.ts + 1
          }
        } else {
          current.ts <- 1
        }

        for (ts in 1:length(trials[[trial.types[j]]])) {
          timesteps <- rbind(timesteps, c(phase.name, current.trial, current.ts, trimws(unlist(strsplit(trials[[trial.types[j]]][ts], ",")))))
          current.ts <- current.ts + 1
        }
        current.trial <- current.trial + 1
        current.ts <- 1
      }
    }
  }

  columns <- suppressWarnings(which(!is.na(sapply(timesteps[1, ], as.numeric)) == T))

  timesteps[, columns] <- as.data.frame(apply(timesteps[, columns], 2, as.numeric))

  colnames(timesteps)[1:3] <- c("Phase", "Trial", "TimeStep")

  return(timesteps)
}

create.NPEs <- function(npe) {
  NPEs <- as.data.frame(matrix(nrow = length(npe), ncol = 9))

  for (n in 1:length(npe)) {
    NPEs[n, ] <- trimws(unlist(strsplit(npe[n], ",")))
  }
  columns <- suppressWarnings(which(!is.na(sapply(NPEs[1, ], as.numeric)) == T))

  NPEs[, columns] <- as.data.frame(apply(NPEs[, columns], 2, as.numeric))
  colnames(NPEs) <- c("NPE", "Type", "Layer", "Activation", "Temporal.Summation", "Activation.Decay", "mu", "sigma", "logisSigma")
  return(NPEs)
}

create.Connections <- function(conn, NPEs) {
  Connections <- as.data.frame(matrix(nrow = length(conn), ncol = 7))

  for (n in 1:length(conn)) {
    Connections[n, ] <- trimws(unlist(strsplit(conn[n], ",")))

    if (any(!Connections[n, 1:2] %in% NPEs)) stop("Either the preSinaptic or the postSinaptic NPE is not a member of the NPEs")
  }
  columns <- suppressWarnings(which(!is.na(sapply(Connections[1, ], as.numeric)) == T))

  Connections[, columns] <- as.data.frame(apply(Connections[, columns], 2, as.numeric))
  colnames(Connections) <- c("PreSinapticNPE", "PostSinapticNPE", "Weight", "alpha", "beta", "alpha_prime", "beta_prime")
  return(Connections)
}

verify.Phases <- function(phases, trials) {
  # phases is a character vector with comma delimited characters. For example,
  # "Training 1, Random,A+/X-,100-100,False,30,30,ITI entrenamiento"
  # "Training 2, In bulk,AX+,100,False,30,30,ITI entrenamiento"
  # "Test,       In bulk,Prueba X,25,False,30,30,ITI entrenamiento"

  # trials is a list with named elements
  # each element is a named vector of comma delimited characters
  # for example, trials[["A+"]] might be:

  #   [1] "S.1,1,S.2,0.65,S.3,0,US,0,True"
  #   [2] "S.1,1,S.2,0.65,S.3,0,US,0,True"
  #   [3] "S.1,1,S.2,0.65,S.3,0,US,0,True"
  #   [4] "S.1,1,S.2,0.65,S.3,0,US,0,True"
  #   [5] "S.1,1,S.2,0.65,S.3,0,US,1,True"

  for (i in 1:length(phases)) {
    current.ts <- 1
    current.trial <- 1

    current.phase <- trimws(unlist(strsplit(phases[i], ",")))

    phase.name <- current.phase[1]
    trial.order <- tolower(current.phase[2])
    trial.types <- trimws(unlist(strsplit(current.phase[3], "/")))

    if (any(!(trial.types %in% names(trials)))) stop("One of more trial names do not match trial names in phases")
  }

  cat("No errors found in the specification of the contingencies")
}

#################################################################################################

ui <- function(request) {
  tagList(
    useShinyjs(),
    tags$script('$(function () { $("[data-toggle=\'tooltip\']").tooltip(); })'),
    tags$head(
      tags$style(HTML("
    /* Styles for the splash screen */
    :root {
        --ddm-bg: #eef3f7;
        --ddm-panel: #ffffff;
        --ddm-accent: #1f6b87;
        --ddm-accent-soft: #e7f2f7;
        --ddm-border: #d7e2eb;
    }
    #splash-screen {
        position: fixed;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background: linear-gradient(135deg, #184e77 0%, #1f7a8c 100%);
        display: flex;
        justify-content: center;
        align-items: center;
        z-index: 9999;
        transition: opacity 0.5s ease-out;
    }
    #splash-content {
        text-align: center;
        color: white;
    }
    #splash-logo {
        width: 200px;
        height: 200px;
        margin-bottom: 30px;
    }
    #splash-title {
        font-size: 3em;
        margin-bottom: 15px;
        opacity: 0;
        transform: translateY(20px);
        transition: opacity 0.5s ease-out, transform 0.5s ease-out;
    }
    #splash-subtitle, .splash-info {
        font-size: 1.2em;
        opacity: 0;
        transform: translateY(20px);
        transition: opacity 0.5s ease-out, transform 0.5s ease-out;
    }
    .node {
        fill: #4CAF50;
    }
    .link {
        stroke: #FFFFFF;
        stroke-width: 2;
    }
    @keyframes fadeIn {
        from { opacity: 0; }
        to { opacity: 1; }
    }
    #splash-logo {
        opacity: 0;
        animation: fadeIn 1s ease-out forwards;
    }
    /* Styles for the main content */
    #main-content {
        display: none;
    }
    .content-wrapper, .right-side {
        background-color: var(--ddm-bg);
        font-family: 'Avenir Next', 'Segoe UI', 'Helvetica Neue', sans-serif;
    }
    .content {
        padding-top: 14px;
    }
    .ux-topbar {
        margin-bottom: 12px;
        position: sticky;
        top: 50px;
        z-index: 900;
    }
    .ux-flow {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
        gap: 8px;
        margin-bottom: 8px;
    }
    .ux-step {
        display: flex;
        align-items: center;
        gap: 8px;
        border: 1px solid var(--ddm-border);
        background: #ffffff;
        border-radius: 10px;
        padding: 8px 10px;
        min-height: 42px;
    }
    .ux-step-num {
        display: inline-flex;
        justify-content: center;
        align-items: center;
        width: 22px;
        height: 22px;
        border-radius: 50%;
        font-size: 12px;
        font-weight: 700;
        border: 1px solid #aac1d1;
        color: #1f4e66;
        background: #f7fbff;
        flex-shrink: 0;
    }
    .ux-step-label {
        font-size: 12px;
        font-weight: 700;
        color: #25455d;
        line-height: 1.2;
    }
    .ux-step.done {
        background: #f0f8f4;
        border-color: #bfe0cc;
    }
    .ux-step.done .ux-step-num {
        background: #d4f0de;
        border-color: #8ec8a5;
        color: #155d3a;
    }
    .ux-step.active {
        border-color: #88b8d1;
        background: #eaf4fb;
        box-shadow: 0 2px 8px rgba(31, 107, 135, 0.15);
    }
    .ux-step.active .ux-step-num {
        background: #1f6b87;
        border-color: #1f6b87;
        color: #ffffff;
    }
    .ux-hint {
        border: 1px solid #cde0ec;
        background: #ffffff;
        border-radius: 10px;
        padding: 10px 12px;
        margin-bottom: 10px;
        color: #1f3f57;
        display: flex;
        align-items: flex-start;
        gap: 10px;
    }
    .ux-hint-icon {
        width: 24px;
        height: 24px;
        border-radius: 50%;
        background: #eaf4fb;
        border: 1px solid #c8dfee;
        color: #1f6b87;
        font-size: 14px;
        font-weight: 700;
        display: inline-flex;
        justify-content: center;
        align-items: center;
        flex-shrink: 0;
    }
    .ux-hint-title {
        font-size: 13px;
        font-weight: 800;
        color: #163850;
        margin-bottom: 2px;
    }
    .ux-hint-text {
        font-size: 13px;
        margin-bottom: 0;
    }
    .box {
        border-top: 0;
        border-radius: 12px;
        box-shadow: 0 8px 20px rgba(25, 50, 75, 0.08);
        border: 1px solid var(--ddm-border);
    }
    .box-title {
        font-weight: 700;
        color: #17324d;
    }
    .box-header {
        border-bottom: 1px solid #edf2f7;
    }
    .sidebar-menu > li > a {
        font-weight: 600;
    }
    .sidebar-menu > li.active > a {
        border-left-color: var(--ddm-accent);
        background: rgba(31, 107, 135, 0.14);
    }
    .btn-primary, .btn-info {
        background-color: var(--ddm-accent);
        border-color: var(--ddm-accent);
    }
    .btn-default {
        border-radius: 8px;
        border: 1px solid #ccd9e5;
        background: #ffffff;
        color: #1f3e56;
        font-weight: 700;
    }
    .btn-default:hover {
        background: #f4f8fb;
        border-color: #afc7d8;
        color: #17364d;
    }
    .box .action-button {
        border-radius: 8px;
        font-weight: 700;
        margin-bottom: 8px;
        max-width: 100%;
        white-space: normal;
    }
    .vis-network .vis-button {
        margin: 0 !important;
        width: 40px !important;
        height: 40px !important;
        min-width: 40px !important;
        min-height: 40px !important;
        border-radius: 50% !important;
        font-size: 18px !important;
        line-height: 38px !important;
    }
    #add_connection {
        width: 100%;
        white-space: normal;
        line-height: 1.2;
        min-height: 44px;
    }
    #reorganize_network {
        display: inline-block;
        width: auto;
        max-width: 100%;
        min-height: 40px;
        padding: 8px 14px;
        box-sizing: border-box;
        white-space: nowrap;
    }
    .btn-block {
        text-align: left;
    }
    .modal-content {
        border-radius: 12px;
        border: 1px solid #d3e2ee;
    }
    .modal-header {
        border-bottom: 1px solid #edf2f7;
    }
    .modal-footer {
        border-top: 1px solid #edf2f7;
    }
    .dataTables_wrapper .dataTables_filter input {
        border-radius: 8px;
        border: 1px solid #c9d7e3;
    }
    .sidebar .radio label {
        color: #d9e6f0;
        font-weight: 600;
    }
    #add_npe, #add_connection, #add_trial, #add_iti, #add_contingency, #run_simulation, #create_architecture,
    #graficar, #graficar_general, #load_extinction_template {
        background-color: #1f6b87;
        border-color: #1f6b87;
        color: #ffffff;
    }
    #add_npe:hover, #add_connection:hover, #add_trial:hover, #add_iti:hover, #add_contingency:hover, #run_simulation:hover, #create_architecture:hover,
    #graficar:hover, #graficar_general:hover, #load_extinction_template:hover {
        background-color: #16566d;
        border-color: #16566d;
        color: #ffffff;
    }
    .form-group > label, .control-label {
        font-size: 13px;
        font-weight: 700;
        color: #1a3a52;
    }
    .form-control, .selectize-input {
        border-radius: 8px;
        border-color: #c9d7e3;
        box-shadow: none;
    }
    .selectize-input.focus, .form-control:focus {
        border-color: #86b6ce;
        box-shadow: 0 0 0 2px rgba(31, 107, 135, 0.12);
    }
    .network-canvas {
        border: 1px solid var(--ddm-border);
        border-radius: 10px;
        background: #ffffff;
        padding: 6px;
    }
    .network-actions {
        margin-top: 10px;
        display: flex;
        gap: 10px;
        flex-wrap: wrap;
    }
    @media (max-width: 900px) {
      .ux-flow {
        grid-template-columns: repeat(2, minmax(120px, 1fr));
      }
    }
    /* Style for help icons */
    .help-icon {
        color: #3c8dbc;
        margin-left: 5px;
        cursor: pointer;
    }
    /* Styles for tables in the help section */
    .table-bordered {
      border: 1px solid #ddd;
    }
    .table-bordered > thead > tr > th,
    .table-bordered > tbody > tr > td {
      border: 1px solid #ddd;
    }
    .table-striped > tbody > tr:nth-of-type(odd) {
      background-color: #f9f9f9;
    }
    .table-hover > tbody > tr:hover {
      background-color: #f5f5f5;
    }
  "))
    ),
    # Splash screen
    div(
      id = "splash-screen",
      div(
        id = "splash-content",
        tags$svg(
          id = "splash-logo", viewBox = "0 0 100 100",
          tags$g(
            tags$circle(class = "node", cx = "50", cy = "20", r = "5"),
            tags$circle(class = "node", cx = "20", cy = "50", r = "5"),
            tags$circle(class = "node", cx = "80", cy = "50", r = "5"),
            tags$circle(class = "node", cx = "35", cy = "80", r = "5"),
            tags$circle(class = "node", cx = "65", cy = "80", r = "5"),
            tags$line(class = "link", x1 = "50", y1 = "20", x2 = "20", y2 = "50"),
            tags$line(class = "link", x1 = "50", y1 = "20", x2 = "80", y2 = "50"),
            tags$line(class = "link", x1 = "20", y1 = "50", x2 = "35", y2 = "80"),
            tags$line(class = "link", x1 = "80", y1 = "50", x2 = "65", y2 = "80"),
            tags$line(class = "link", x1 = "35", y1 = "80", x2 = "65", y2 = "80")
          )
        ),
        h1(id = "splash-title", "Diffuse Discrepancy Model"),
        p(id = "splash-subtitle", "Based on Donahoe, Burgos and Palmer (1993)"),
        p(class = "splash-info", "Interface designed by Miguel Ángel Aguayo Mendoza"),
        p(class = "splash-info", "Last modified: February 2026"),
        p(class = "splash-info", "University of Guadalajara")
      )
    ),
    # Main content of the application
    div(
      id = "main-content",
      dashboardPage(
        dashboardHeader(title = "DiffDiscM Simulator"),
        dashboardSidebar(
          sidebarMenu(
            id = "main_nav",
            menuItem("Home", tabName = "home", icon = icon("home")),
            menuItem("Network Architecture", tabName = "network", icon = icon("project-diagram")),
            menuItem("Create Trials", tabName = "trials", icon = icon("list")),
            menuItem("Configure Contingencies", tabName = "contingencies", icon = icon("cogs")),
            menuItem("Simulate", tabName = "simulate", icon = icon("play")),
            menuItem("Individual Results", tabName = "results_individual", icon = icon("chart-line")),
            menuItem("General Results", tabName = "results_general", icon = icon("chart-bar")),
            menuItem("Help", tabName = "help", icon = icon("question-circle"))
          ),
          tags$div(
            style = "padding: 0 10px 10px 10px;",
            radioButtons(
              "app_mode",
              "Mode",
              choices = c("Beginner", "Advanced"),
              selected = "Beginner"
            )
          ),
          tags$div(
            style = "position: absolute; bottom: 0; left: 0; right: 0.4; padding: 10px;",
            actionButton("close_app", "Close Program",
              icon = icon("power-off"),
              style = "width: 100%; color: #fff; background-color: #d9534f; border-color: #d43f3a;"
            )
          )
        ),
        dashboardBody(
          div(
            class = "ux-topbar",
            uiOutput("workflow_progress"),
            uiOutput("screen_hint")
          ),
          tabItems(
            # Home tab
            tabItem(
              tabName = "home",
              fluidRow(
                box(
                  title = "Welcome to the Diffuse Discrepancy Model (DiffDiscM) Simulator",
                  width = 12,
                  p("The Diffuse Discrepancy Model (DiffDiscM) is a powerful tool for simulating learning and conditioning phenomena in behavioral sciences."),
                  p("Originally developed by Donahoe, Burgos and Palmer (1993), the DiffDiscM offers a connectionist interpretation of the unified principle of reinforcement for operant and Pavlovian conditioning."),
                  h4("Key features:"),
                  tags$ul(
                    tags$li("Simulates both Pavlovian and operant conditioning."),
                    tags$li("Based on principles of neuroanatomy and neurophysiology."),
                    tags$li("Uses activation and learning rules to model the behavior of neural processing units (NPUs)."),
                    tags$li("Incorporates hippocampal and dopaminergic systems in the learning process.")
                  ),
                  h4("How to use this simulator:"),
                  tags$ol(
                    tags$li("Set up the network architecture in the 'Network Architecture' tab."),
                    tags$li("Define trials in the 'Create Trials' tab."),
                    tags$li("Configure contingencies in the 'Configure Contingencies' tab."),
                    tags$li("Run the simulation in the 'Simulate' tab."),
                    tags$li("Analyze results in the 'Individual Results' and 'General Results' sections.")
                  ),
                  p("For more information on how to use each component of the simulator, please refer to the 'Help' section."),
                  hr(),
                  p("Quick start for new users: load a minimal extinction template (architecture + trials + contingencies)."),
                  actionButton("load_extinction_template", "Load Minimal Extinction Template", icon = icon("flask")),
                  br(), br(),
                  p(strong("Storage shortcuts (no explorer required):")),
                  actionButton("use_workdir_paths", "Use Current Folder In All Path Fields", icon = icon("folder-open")),
                  actionButton("use_downloads_paths", "Use Downloads In All Path Fields", icon = icon("download")),
                  br(), br(),
                  textOutput("working_dir_label"),
                  textOutput("downloads_dir_label")
                )
              )
            ),
            # Network architecture tab
            tabItem(
              tabName = "network",
              fluidRow(
                box(
                  title = "Units (NPUs)",
                  width = 6,
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("npe_name", "NPU Name"),
                    icon("question-circle", id = "help_npe_name", class = "help-icon")
                  ),
                  fluidRow(
                    column(
                      6,
                      div(
                        style = "display: flex; align-items: center;",
                        selectInput("npe_type", "NPU Type",
                          choices = c("Excitatory", "Inhibitory"),
                          width = "100%"
                        ),
                        icon("question-circle", id = "help_npe_type", class = "help-icon")
                      )
                    ),
                    column(
                      6,
                      div(
                        style = "display: flex; align-items: center;",
                        selectInput("npe_layer", "Layer",
                          choices = c("PrimarySensory", "AssociativeSensory", "Hippocampal", "AssociativeMotor", "PrimaryMotor", "Dopaminergic"),
                          width = "100%"
                        ),
                        icon("question-circle", id = "help_npe_layer", class = "help-icon")
                      )
                    )
                  ),
                  fluidRow(
                    column(4, numericInput("npe_activation", "Initial Activation", value = 0, min = 0, max = 1, step = 0.1)),
                    column(4, numericInput("npe_temporal_summation", "Temporal Summation (τ)", value = 0.1, min = 0, max = 1, step = 0.1)),
                    column(4, numericInput("npe_decay", "Decay Rate (κ)", value = 0.1, min = 0, max = 1, step = 0.1))
                  ),
                  div(
                    id = "advanced_npu_params",
                    fluidRow(
                      column(4, numericInput("npe_mu", "Threshold Mean (μ)", value = 0.2, min = 0, max = 1, step = 0.1)),
                      column(4, numericInput("npe_sigma", "Threshold Deviation (σ)", value = 0.15, min = 0, max = 1, step = 0.1)),
                      column(4, numericInput("npe_logistic_slope", "Logistic Slope", value = 0.1, min = 0, max = 1, step = 0.1))
                    )
                  ),
                  actionButton("add_npe", "Add NPU", icon = icon("plus-circle")),
                  hr(),
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("npe_file_name", "NPUs File Name"),
                    icon("question-circle", id = "help_npe_file_name", class = "help-icon")
                  ),
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("npe_file_path", "NPUs Directory Path (optional)"),
                    icon("question-circle", id = "help_npe_file_path", class = "help-icon")
                  ),
                  actionButton("save_npes", "Save NPUs", icon = icon("save")),
                  actionButton("import_npes", "Import NPUs", icon = icon("folder-open"))
                ),
                box(
                  title = "Connections",
                  width = 6,
                  div(
                    style = "display: flex; align-items: center;",
                    selectInput("conn_pre", "Source NPU", choices = NULL),
                    icon("question-circle", id = "help_conn_pre", class = "help-icon")
                  ),
                  div(
                    style = "display: flex; align-items: center;",
                    selectInput("conn_post", "Target NPU", choices = NULL),
                    icon("question-circle", id = "help_conn_post", class = "help-icon")
                  ),
                  fluidRow(
                    column(4, numericInput("conn_weight", "Initial Weight", value = 0.1, min = 0, max = 1, step = 0.1)),
                    column(4, div(id = "conn_alpha_wrap", numericInput("conn_alpha", "Increment Rate (α)", value = 0.5, min = 0, max = 1, step = 0.1))),
                    column(4, div(id = "conn_beta_wrap", numericInput("conn_beta", "Decrement Rate (β)", value = 0.12, min = 0, max = 1, step = 0.1)))
                  ),
                  fluidRow(
                    column(6, div(id = "conn_alpha_prime_wrap", numericInput("conn_alpha_prime", "Inhibitory Increment Rate (α')", value = 0.5, min = 0, max = 1, step = 0.1))),
                    column(6, div(id = "conn_beta_prime_wrap", numericInput("conn_beta_prime", "Inhibitory Decrement Rate (β')", value = 0.12, min = 0, max = 1, step = 0.1)))
                  ),
                  fluidRow(
                    column(
                      6,
                      actionButton("add_connection", "Add Connection", icon = icon("link"), style = "width: 100%;"),
                      icon("question-circle", id = "help_add_connection", class = "help-icon")
                    ),
                    column(
                      6,
                      div(
                        id = "reminder_us_d",
                        style = "background-color: #fff3cd; border: 1px solid #ffeeba; border-radius: 4px; padding: 6px 12px; cursor: help;",
                        span(
                          "Reminder",
                          style = "font-weight: bold; color: #856404;",
                          `data-toggle` = "tooltip",
                          `data-placement` = "top",
                          title = "The connection between US and D must have an exact maximum weight of 1. Don't forget to press 'Create Architecture' at the end of this section."
                        )
                      )
                    )
                  ),
                  hr(),
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("conn_file_name", "Connections File Name"),
                    icon("question-circle", id = "help_conn_file_name", class = "help-icon")
                  ),
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("conn_file_path", "Connections Directory Path (optional)"),
                    icon("question-circle", id = "help_conn_file_path", class = "help-icon")
                  ),
                  actionButton("save_connections", "Save Connections", icon = icon("save")),
                  actionButton("import_connections", "Import Connections", icon = icon("folder-open"))
                )
              ),
              fluidRow(
                box(
                  title = "Network Visualization",
                  width = 12,
                  div(
                    class = "network-canvas",
                    visNetworkOutput("network_plot", height = "520px")
                  ),
                  div(
                    class = "network-actions",
                    actionButton("reorganize_network", "Reset Layout", icon = icon("project-diagram"))
                  )
                )
              ),
              fluidRow(
                box(
                  title = "Defined Units",
                  width = 6,
                  DTOutput("npe_table")
                ),
                box(
                  title = "Defined Connections",
                  width = 6,
                  DTOutput("connection_table")
                )
              ),
              fluidRow(
                box(
                  title = "Create Architecture",
                  width = 12,
                  div(
                    style = "display: flex; align-items: center;",
                    actionButton("create_architecture", "Create Architecture", icon = icon("project-diagram"))
                  )
                )
              )
            ),
            tabItem(
              tabName = "trials",
              fluidRow(
                box(
                  title = "Trial Configuration",
                  width = 12,
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("trial_name", "Trial Type Name"),
                    icon("question-circle", id = "help_trial_name", class = "help-icon")
                  ),
                  div(
                    style = "display: flex; align-items: center;",
                    numericInput("num_moments", "Number of Time Moments",
                      value = 1, min = 1, max = 10, step = 1
                    ),
                    icon("question-circle", id = "help_num_moments", class = "help-icon")
                  ),
                  uiOutput("dynamic_inputs"),
                  checkboxInput("learning_rule", "Active Learning Rule", value = TRUE),
                  actionButton("add_trial", "Add Trial", icon = icon("plus-circle")),
                  actionButton("add_iti", "Add ITI", icon = icon("clock"))
                )
              ),
              fluidRow(
                box(
                  title = "Created Trials",
                  width = 12,
                  uiOutput("trial_buttons")
                )
              ),
              fluidRow(
                box(
                  title = "Save/Import Trials",
                  width = 12,
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("trials_file_name", "Trials File Name (without extension)"),
                    icon("question-circle", id = "help_trials_file_name", class = "help-icon")
                  ),
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("trials_file_path", "Trials Directory Path (optional)"),
                    icon("question-circle", id = "help_trials_file_path", class = "help-icon")
                  ),
                  actionButton("save_trials", "Save Trials", icon = icon("save")),
                  actionButton("import_trials", "Import Trials", icon = icon("folder-open"))
                )
              )
            ),
            tabItem(
              tabName = "contingencies",
              fluidRow(
                box(
                  title = "Configure Contingencies",
                  width = 12,
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("phase_name", "Phase or Condition Name"),
                    icon("question-circle", id = "help_phase_name", class = "help-icon")
                  ),
                  div(
                    style = "display: flex; align-items: center;",
                    selectInput("presentation_mode", "Trial Presentation Mode",
                      choices = c("Random", "In bulk", "Alternated")
                    ),
                    icon("question-circle", id = "help_presentation_mode", class = "help-icon")
                  ),
                  div(
                    style = "display: flex; align-items: center;",
                    selectInput("trial_types", "Trial Types to Present",
                      choices = NULL, multiple = TRUE
                    ),
                    icon("question-circle", id = "help_trial_types", class = "help-icon")
                  ),
                  uiOutput("trial_numbers"),
                  checkboxInput("reset_activations", "Reset Activations", value = TRUE),
                  conditionalPanel(
                    condition = "!input.reset_activations",
                    numericInput("min_iti", "Minimum ITI Value", value = 30, min = 1),
                    numericInput("max_iti", "Maximum ITI Value", value = 30, min = 1),
                    div(
                      style = "display: flex; align-items: center;",
                      selectInput("iti_trial", "Add Created ITI", choices = NULL),
                      icon("question-circle", id = "help_iti_trial", class = "help-icon")
                    )
                  ),
                  actionButton("add_contingency", "Add Contingency", icon = icon("plus-circle"))
                )
              ),
              fluidRow(
                box(
                  title = "Created Contingencies",
                  width = 12,
                  DTOutput("contingencies_table")
                )
              ),
              fluidRow(
                box(
                  title = "Save/Import Contingencies",
                  width = 12,
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("contingencies_file_name", "Contingencies File Name"),
                    icon("question-circle", id = "help_contingencies_file_name", class = "help-icon")
                  ),
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("contingencies_file_path", "Contingencies Directory Path (optional)"),
                    icon("question-circle", id = "help_contingencies_file_path", class = "help-icon")
                  ),
                  actionButton("save_contingencies", "Save Contingencies", icon = icon("save")),
                  actionButton("import_contingencies", "Import Contingencies", icon = icon("folder-open"))
                )
              )
            ),
            tabItem(
              tabName = "simulate",
              fluidRow(
                box(
                  title = "Run Simulation",
                  width = 12,
                  div(
                    style = "display: flex; align-items: center;",
                    numericInput("num_simulations", "Number of Networks", value = 1, min = 1, step = 1),
                    icon("question-circle", id = "help_num_simulations", class = "help-icon")
                  ),
                  actionButton("run_simulation", "Run Simulation", icon = icon("play")),
                  actionButton("save_networks", "Save Networks", icon = icon("download")),
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("sim_file_name", "Simulation File Name"),
                    icon("question-circle", id = "help_sim_file_name", class = "help-icon")
                  ),
                  div(
                    style = "display: flex; align-items: center;",
                    textInput("sim_file_path", "Simulation Directory Path (optional)"),
                    icon("question-circle", id = "help_sim_file_path", class = "help-icon")
                  ),
                  actionButton("save_simulation", "Save Simulation", icon = icon("save")),
                  actionButton("load_simulation", "Load Simulation", icon = icon("upload")),
                  progressBar("sim_progress", value = 0, display_pct = TRUE),
                  verbatimTextOutput("simulation_status")
                )
              ),
              fluidRow(
                box(
                  title = "Simulation Results",
                  width = 12,
                  DTOutput("simulation_results")
                )
              )
            ),
            tabItem(
              tabName = "results_individual",
              fluidRow(
                box(
                  title = "Select Simulation",
                  width = 12,
                  selectInput("selected_simulation", "Select Network", choices = NULL)
                )
              ),
              fluidRow(
                box(
                  title = "Plot Results",
                  width = 12,
                  actionButton("graficar", "Plot Results", icon = icon("line-chart"))
                )
              ),
              fluidRow(
                box(
                  title = "Activation Plot",
                  width = 12,
                  selectizeInput("activations_units", "Select Units",
                    choices = NULL, multiple = TRUE
                  ),
                  selectInput("selected_timestep", "Select Time Step", choices = NULL),
                  plotlyOutput("activations_plot", height = "400px")
                )
              ),
              fluidRow(
                box(
                  title = "Connection Weights Plot",
                  width = 12,
                  selectizeInput("weights_connections", "Select Connections",
                    choices = NULL, multiple = TRUE
                  ),
                  plotlyOutput("weights_plot", height = "400px")
                )
              ),
              fluidRow(
                box(
                  title = "Aggregate Measures",
                  width = 12,
                  selectInput("aggregate_phase", "Select Phase", choices = NULL),
                  selectizeInput("aggregate_units", "Select Units", choices = NULL, multiple = TRUE),
                  selectInput("aggregate_timestep", "Select Time Step", choices = NULL),
                  selectInput("aggregate_measure", "Select Measure",
                    choices = c("Mean" = "mean", "Median" = "median")
                  ),
                  selectInput("aggregate_error", "Select Error",
                    choices = c("Standard Error" = "se", "Standard Deviation" = "sd")
                  ),
                  plotlyOutput("aggregate_plot", height = "400px")
                )
              )
            ),
            tabItem(
              tabName = "results_general",
              fluidRow(
                box(
                  title = "Plot General Results",
                  width = 12,
                  actionButton("graficar_general", "Plot General Results", icon = icon("bar-chart"))
                )
              ),
              fluidRow(
                box(
                  title = "Activation Plot (All Networks)",
                  width = 12,
                  selectizeInput("activations_units_general", "Select Units",
                    choices = NULL, multiple = TRUE
                  ),
                  selectInput("selected_timestep_general", "Select Time Step", choices = NULL),
                  plotlyOutput("activations_plot_general", height = "600px")
                )
              ),
              fluidRow(
                box(
                  title = "Connection Weights Plot (All Networks)",
                  width = 12,
                  selectizeInput("weights_connections_general", "Select Connections",
                    choices = NULL, multiple = TRUE
                  ),
                  plotlyOutput("weights_plot_general", height = "600px")
                )
              ),
              fluidRow(
                box(
                  title = "Aggregate Measures (All Networks)",
                  width = 12,
                  selectInput("aggregate_phase_general", "Select Phase", choices = NULL),
                  selectizeInput("aggregate_units_general", "Select Units", choices = NULL, multiple = TRUE),
                  selectInput("aggregate_timestep_general", "Select Time Step", choices = NULL),
                  selectInput("aggregate_measure_general", "Select Measure",
                    choices = c("Mean" = "mean", "Median" = "median")
                  ),
                  selectInput("aggregate_error_general", "Select Error",
                    choices = c("Standard Error" = "se", "Standard Deviation" = "sd")
                  ),
                  plotlyOutput("aggregate_plot_general", height = "400px")
                )
              )
            ),
            tabItem(
              tabName = "help",
              fluidRow(
                box(
                  title = "User Guide",
                  width = 12,
                  tabsetPanel(
                    tabPanel(
                      "Introduction",
                      h4("Welcome to the DiffDiscM Simulator"),
                      p("This guide will help you effectively use the Diffuse Discrepancy Model (DiffDiscM) Simulator."),
                      h4("Model features:"),
                      tags$ul(
                        tags$li("Based on principles of neuroanatomy and neurophysiology."),
                        tags$li("Uses activation and learning rules to model the behavior of neural processing units (NPUs)."),
                        tags$li("Incorporates hippocampal and dopaminergic systems in the learning process."),
                        tags$li("Allows simulation of both Pavlovian and operant conditioning.")
                      )
                    ),
                    tabPanel(
                      "Network Architecture",
                      h4("Neural Network Configuration"),
                      p("In this section, you can define the structure of your neural network:"),
                      tags$ul(
                        tags$li("Add units (NPUs) specifying their characteristics."),
                        tags$li("Establish connections between units."),
                        tags$li("Visualize the network to verify its structure.")
                      ),
                      h4("NPU Types:"),
                      tags$ul(
                        tags$li("PrimarySensory: Simulate primary sensory effects of environmental events."),
                        tags$li("AssociativeSensory: Represent associative sensory areas."),
                        tags$li("Hippocampal: Simulate hippocampal areas involved in conditioning."),
                        tags$li("AssociativeMotor: Represent associative motor areas."),
                        tags$li("PrimaryMotor: Simulate primary motor precursors."),
                        tags$li("Dopaminergic: Simulate dopaminergic areas like the ventral tegmental area."),
                        tags$li("US (R*): Simulates the unconditioned stimulus or reinforcer. This unit is crucial for learning and conditioning in the model.")
                      ),
                      h4("US (R*) Unit:"),
                      p("The US (Unconditioned Stimulus) or R* unit is a special unit in the model that simulates the unconditioned stimulus or reinforcer. Some important characteristics of this unit are:"),
                      tags$ul(
                        tags$li("Represents biologically significant events, such as food in conditioning experiments."),
                        tags$li("Has a fixed connection with maximum weight to the dopaminergic (D) unit. Therefore, you must set the weight of the US to D connection to 1."),
                        tags$li("Its activation produces an unconditioned response and also serves as a reinforcement signal for learning."),
                        tags$li("It is fundamental for simulating both Pavlovian and operant conditioning.")
                      ),
                      p("The D and US units are already created by default, just for you to connect them. At the moment, it is not possible to add more than 1 US unit."),
                      h4("Example of NPUs file (CSV):"),
                      div(
                        style = "overflow-x: auto;",
                        tags$table(
                          class = "table table-bordered table-striped table-hover",
                          tags$thead(
                            tags$tr(
                              tags$th("NPU"), tags$th("Type"), tags$th("Layer"),
                              tags$th("Activation"), tags$th("Temporal.Summation"),
                              tags$th("Activation.Decay"), tags$th("mu"),
                              tags$th("sigma"), tags$th("logisSigma"),
                              tags$th("x"), tags$th("y")
                            )
                          ),
                          tags$tbody(
                            tags$tr(
                              tags$td("US"), tags$td("Excitatory"), tags$td("US"),
                              tags$td("0"), tags$td("0.1"), tags$td("0.1"),
                              tags$td("0.2"), tags$td("0.15"), tags$td("0.1"),
                              tags$td("0.86647756"), tags$td("0.77925187")
                            ),
                            tags$tr(
                              tags$td("D"), tags$td("Excitatory"), tags$td("Dopaminergic"),
                              tags$td("0"), tags$td("0.1"), tags$td("0.1"),
                              tags$td("0.2"), tags$td("0.15"), tags$td("0.1"),
                              tags$td("0.4778924"), tags$td("0.72927362")
                            ),
                            tags$tr(
                              tags$td("S.1"), tags$td("Excitatory"), tags$td("PrimarySensory"),
                              tags$td("0"), tags$td("0.1"), tags$td("0.1"),
                              tags$td("0.2"), tags$td("0.15"), tags$td("0.1"),
                              tags$td("0.72974392"), tags$td("0.94505817")
                            )
                          )
                        )
                      ),
                      h4("Example of Connections file (CSV):"),
                      div(
                        style = "overflow-x: auto;",
                        tags$table(
                          class = "table table-bordered table-striped table-hover",
                          tags$thead(
                            tags$tr(
                              tags$th("PreSynapticNPU"), tags$th("PostSynapticNPU"),
                              tags$th("Weight"), tags$th("alpha"), tags$th("beta"),
                              tags$th("alpha_prime"), tags$th("beta_prime")
                            )
                          ),
                          tags$tbody(
                            tags$tr(
                              tags$td("US"), tags$td("D"), tags$td("1"),
                              tags$td("0.5"), tags$td("0.12"), tags$td("0.5"), tags$td("0.12")
                            ),
                            tags$tr(
                              tags$td("S.1"), tags$td("S.2"), tags$td("0.1"),
                              tags$td("0.5"), tags$td("0.12"), tags$td("0.5"), tags$td("0.12")
                            ),
                            tags$tr(
                              tags$td("S.2"), tags$td("H1"), tags$td("0.1"),
                              tags$td("0.5"), tags$td("0.12"), tags$td("0.5"), tags$td("0.12")
                            )
                          )
                        )
                      ),
                      p("Note: The NPU names are examples. You can use any name you want for your units.")
                    ),
                    tabPanel(
                      "Create Trials",
                      h4("Trial Design"),
                      p("Here you can configure the different types of trials:"),
                      tags$ul(
                        tags$li("Define trial types and their characteristics."),
                        tags$li("Specify the time moments for each trial."),
                        tags$li("Create Inter-Trial Interval (ITI) trials if necessary.")
                      ),
                      p("Each trial is defined as a series of NPU activations at different time moments."),
                      p("The trials are saved in a file with .rds format. This is an R binary format and cannot be directly viewed as text. It contains a list structure with the trials you have defined in the interface.")
                    ),
                    tabPanel(
                      "Contingency Configuration",
                      h4("Contingency Design"),
                      p("In this section, you can configure the experimental contingencies:"),
                      tags$ul(
                        tags$li("Define experimental phases."),
                        tags$li("Specify the trial presentation mode (Random, In bulk, Alternated). Legacy files using 'In block' are also accepted."),
                        tags$li("Configure inter-trial intervals (ITI) if necessary.")
                      ),
                      p("The contingencies determine how different types of trials are presented during the simulation."),
                      p("Typically, Random is used for training trials and In bulk for test trials."),
                      h4("Example of Contingencies file (CSV):"),
                      div(
                        style = "overflow-x: auto;",
                        tags$table(
                          class = "table table-bordered table-striped table-hover",
                          tags$thead(
                            tags$tr(
                              tags$th("Contingency")
                            )
                          ),
                          tags$tbody(
                            tags$tr(
                              tags$td("Training, Random, S.1/S.2, 100-100, False")
                            ),
                            tags$tr(
                              tags$td("Test, In bulk, TestX, 25, False")
                            )
                          )
                        )
                      ),
                      p("Note: Each line, including the 'Contingency' header, is in quotes in the actual CSV file. The phase names and trial types are examples and can be customized according to your needs.")
                    ),
                    tabPanel(
                      "Simulation",
                      h4("Running Simulations"),
                      p("In this section, you can run your simulations:"),
                      tags$ul(
                        tags$li("Specify the number of networks to simulate."),
                        tags$li("Start the simulation and monitor its progress."),
                        tags$li("Save and load simulations for later analysis.")
                      ),
                      p("During the simulation, the model updates the NPU activations and connection weights according to the activation and learning rules.")
                    ),
                    tabPanel(
                      "Results Analysis",
                      h4("Data Visualization and Analysis"),
                      p("Here you can explore and analyze the results of your simulations:"),
                      tags$ul(
                        tags$li("View activation and connection weight graphs."),
                        tags$li("Compare results between different simulations."),
                        tags$li("Perform statistical analysis of the data.")
                      ),
                      p("The results will allow you to understand how the model simulates different learning and conditioning phenomena.")
                    )
                  )
                )
              ),
              fluidRow(
                box(
                  title = "Technical Support",
                  width = 12,
                  p("If you encounter any problems or have additional questions, please contact:"),
                  tags$ul(
                    tags$li("Miguel Ángel Aguayo Mendoza"),
                    tags$li("Email: miguel.aguayo@academicos.udg.mx"),
                    tags$li("University of Guadalajara")
                  ),
                  p(
                    "For more information, visit the experimental and theoretical research laboratory in Learning, Conditioning and Adaptive Behavior: ",
                    a("Website", href = "http://www.ceic.cucba.udg.mx/Investigacion/laboratorios?id=13", target = "_blank")
                  )
                )
              )
            )
          )
        )
      )
    ),
    # Script to control the splash screen animation
    tags$script(HTML("
      $(document).ready(function() {
        setTimeout(function() {
          $('#splash-title').css({'opacity': '1', 'transform': 'translateY(0)'});
        }, 500);
        setTimeout(function() {
          $('#splash-subtitle').css({'opacity': '1', 'transform': 'translateY(0)'});
        }, 1000);
        setTimeout(function() {
          $('.splash-info').css({'opacity': '1', 'transform': 'translateY(0)'});
        }, 1500);
        setTimeout(function() {
          $('#splash-screen').css('opacity', '0');
        }, 3000);
        setTimeout(function() {
          $('#splash-screen').hide();
          $('#main-content').show();
        }, 3500);
      });
    ")),
    # Bridge clicks on help icons to Shiny inputs
    tags$script(HTML("
      $(document).on('click', '.help-icon', function(e) {
        e.preventDefault();
        var id = $(this).attr('id');
        if (id) {
          Shiny.setInputValue(id, Date.now(), {priority: 'event'});
        }
      });
    ")),
    # Tooltips
    bsTooltip("help_npe_name", "Name the network unit. It is recommended NOT to use apostrophes or quotes, and to use a short name. For example: S.1", placement = "right", trigger = "hover"),
    bsTooltip("help_npe_type", "Select whether the unit is excitatory or inhibitory", placement = "right", trigger = "hover"),
    bsTooltip("help_npe_layer", "Select the layer to which the unit belongs. The dopaminergic unit and US already exist by default", placement = "right", trigger = "hover"),
    bsTooltip("help_npe_file_name", "Name of the file where the NPUs will be saved", placement = "right", trigger = "hover"),
    bsTooltip("help_npe_file_path", "Directory where NPUs are saved/loaded. Optional: leave blank to use the current folder, or use Home shortcuts.", placement = "right", trigger = "hover"),
    bsTooltip("help_conn_pre", "Select the source NPU of the connection", placement = "right", trigger = "hover"),
    bsTooltip("help_conn_post", "Select the target NPU of the connection", placement = "right", trigger = "hover"),
    bsTooltip("help_add_connection", "Add the defined connection to the network", placement = "right", trigger = "hover"),
    bsTooltip("help_conn_file_name", "Name of the file where the connections will be saved", placement = "right", trigger = "hover"),
    bsTooltip("help_conn_file_path", "Directory where connections are saved/loaded. Optional: leave blank to use the current folder, or use Home shortcuts.", placement = "right", trigger = "hover"),
    bsTooltip("help_trial_name", "Identifying name for the trial type", placement = "right", trigger = "hover"),
    bsTooltip("help_num_moments", "Number of time moments in the trial", placement = "right", trigger = "hover"),
    bsTooltip("help_trials_file_name", "Name of the file where the trials will be saved", placement = "right", trigger = "hover"),
    bsTooltip("help_trials_file_path", "Directory where trials are saved/loaded. Optional: leave blank to use the current folder, or use Home shortcuts.", placement = "right", trigger = "hover"),
    bsTooltip("help_phase_name", "Name of the experimental phase or condition", placement = "right", trigger = "hover"),
    bsTooltip("help_presentation_mode", "Mode of presentation of the trials", placement = "right", trigger = "hover"),
    bsTooltip("help_trial_types", "Select the types of trials to present in this phase", placement = "right", trigger = "hover"),
    bsTooltip("help_iti_trial", "Select the ITI trial to add", placement = "right", trigger = "hover"),
    bsTooltip("help_contingencies_file_name", "Name of the file where the contingencies will be saved", placement = "right", trigger = "hover"),
    bsTooltip("help_contingencies_file_path", "Directory where contingencies are saved/loaded. Optional: leave blank to use the current folder, or use Home shortcuts.", placement = "right", trigger = "hover"),
    bsTooltip("help_num_simulations", "Number of networks to simulate", placement = "right", trigger = "hover"),
    bsTooltip("help_sim_file_name", "Name of the file where the simulation will be saved", placement = "right", trigger = "hover"),
    bsTooltip("help_sim_file_path", "Directory where simulations are saved/loaded. Optional: leave blank to use the current folder, or use Home shortcuts.", placement = "right", trigger = "hover")
  )
}
#################################################################################################
# Servidor completo
server <- function(input, output, session) {
  # Inicialización de valores reactivos
  npes <- reactiveVal(data.frame(
    NPE = c("US", "D"),
    Type = c("Excitatory", "Excitatory"),
    Layer = c("US", "Dopaminergic"),
    Activation = c(0, 0),
    Temporal.Summation = c(0.1, 0.1),
    Activation.Decay = c(0.1, 0.1),
    mu = c(0.2, 0.2),
    sigma = c(0.15, 0.15),
    logisSigma = c(0.1, 0.1),
    x = c(-180, 520),
    y = c(320, 270),
    stringsAsFactors = FALSE
  ))

  connections <- reactiveVal(data.frame(
    PreSinapticNPE = character(),
    PostSinapticNPE = character(),
    Weight = numeric(),
    alpha = numeric(),
    beta = numeric(),
    alpha_prime = numeric(),
    beta_prime = numeric(),
    stringsAsFactors = FALSE
  ))

  trials <- reactiveVal(list())
  contingencies <- reactiveVal(character())
  has_iti <- reactiveVal(logical())
  time_steps <- reactiveVal(NULL)
  simulation_results <- reactiveVal(NULL)
  app_mode <- reactiveVal("Beginner")

  path_input_ids <- c(
    "npe_file_path",
    "conn_file_path",
    "trials_file_path",
    "contingencies_file_path",
    "sim_file_path"
  )

  current_workdir <- function() {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }

  downloads_dir <- function() {
    normalizePath(path.expand("~/Downloads"), winslash = "/", mustWork = FALSE)
  }

  resolve_directory <- function(path_value) {
    value <- if (is.null(path_value)) "" else trimws(path_value)
    if (value == "") {
      return(current_workdir())
    }
    normalizePath(path.expand(value), winslash = "/", mustWork = FALSE)
  }

  set_all_path_inputs <- function(path_value) {
    normalized <- normalizePath(path.expand(path_value), winslash = "/", mustWork = FALSE)
    for (id in path_input_ids) {
      updateTextInput(session, id, value = normalized)
    }
  }

  apply_mode <- function(mode) {
    if (mode == "Beginner") {
      shinyjs::hide("advanced_npu_params")
      shinyjs::hide("conn_alpha_wrap")
      shinyjs::hide("conn_beta_wrap")
      shinyjs::hide("conn_alpha_prime_wrap")
      shinyjs::hide("conn_beta_prime_wrap")
      updateNumericInput(session, "npe_mu", value = 0.2)
      updateNumericInput(session, "npe_sigma", value = 0.15)
      updateNumericInput(session, "npe_logistic_slope", value = 0.1)
      updateNumericInput(session, "conn_alpha", value = 0.5)
      updateNumericInput(session, "conn_beta", value = 0.12)
      updateNumericInput(session, "conn_alpha_prime", value = 0.5)
      updateNumericInput(session, "conn_beta_prime", value = 0.12)
    } else {
      shinyjs::show("advanced_npu_params")
      shinyjs::show("conn_alpha_wrap")
      shinyjs::show("conn_beta_wrap")
      shinyjs::show("conn_alpha_prime_wrap")
      shinyjs::show("conn_beta_prime_wrap")
    }
  }

  output$working_dir_label <- renderText({
    paste("Current folder:", current_workdir())
  })

  output$downloads_dir_label <- renderText({
    paste("Downloads:", downloads_dir())
  })

  get_selected_network_index <- function(results_list = NULL) {
    if (is.null(results_list)) {
      results_list <- simulation_results()
    }
    if (is.null(results_list) || length(results_list) == 0) {
      return(NULL)
    }

    selected_label <- input$selected_simulation
    if (is.null(selected_label) || !nzchar(selected_label)) {
      return(1L)
    }

    idx <- suppressWarnings(as.integer(gsub("^Network\\s*", "", selected_label)))
    if (is.na(idx) || idx < 1 || idx > length(results_list)) {
      return(1L)
    }
    idx
  }

  normalize_label <- function(x) {
    if (is.null(x)) {
      return("")
    }
    x <- trimws(as.character(x))
    ascii <- suppressWarnings(iconv(x, from = "", to = "ASCII//TRANSLIT"))
    ascii[is.na(ascii)] <- x[is.na(ascii)]
    tolower(ascii)
  }

  canonical_presentation_mode <- function(mode_value) {
    mode_norm <- gsub("\\s+", " ", normalize_label(mode_value))
    if (mode_norm %in% c("in block", "in blocks", "block", "blocks", "inblock")) {
      return("In bulk")
    }
    if (mode_norm %in% c("in bulk", "inbulk", "bulk")) {
      return("In bulk")
    }
    if (mode_norm %in% c("random")) {
      return("Random")
    }
    if (mode_norm %in% c("alternated", "alternate", "alternating")) {
      return("Alternated")
    }
    trimws(mode_value)
  }

  map_trial_name <- function(label, trial_names) {
    if (is.null(label) || !nzchar(trimws(label))) {
      return(label)
    }
    cleaned <- trimws(label)
    if (cleaned %in% trial_names) {
      return(cleaned)
    }

    normalized_trials <- normalize_label(trial_names)
    normalized_label <- normalize_label(cleaned)
    idx <- which(normalized_trials == normalized_label)

    if (length(idx) == 1) {
      return(trial_names[idx])
    }
    cleaned
  }

  harmonize_contingencies <- function(contingency_vec, trial_list) {
    trial_names <- names(trial_list)
    if (length(trial_names) == 0) {
      stop("No trials were found. Please create or import trials before running simulations.")
    }

    vapply(seq_along(contingency_vec), function(i) {
      raw_parts <- trimws(unlist(strsplit(contingency_vec[i], ",")))
      if (length(raw_parts) < 5) {
        stop(paste("Contingency", i, "is incomplete. It must have at least 5 comma-separated fields."))
      }

      raw_parts[2] <- canonical_presentation_mode(raw_parts[2])

      trial_tokens <- trimws(unlist(strsplit(raw_parts[3], "/")))
      mapped_trials <- vapply(trial_tokens, function(tt) map_trial_name(tt, trial_names), character(1))
      if (any(!(mapped_trials %in% trial_names))) {
        missing_names <- unique(mapped_trials[!(mapped_trials %in% trial_names)])
        stop(
          paste(
            "Contingency", i, "references unknown trial type(s):",
            paste(missing_names, collapse = ", ")
          )
        )
      }
      raw_parts[3] <- paste(mapped_trials, collapse = "/")

      has_iti <- as.logical(trimws(raw_parts[5]))
      if (!is.na(has_iti) && has_iti && length(raw_parts) >= 8) {
        mapped_iti <- map_trial_name(raw_parts[8], trial_names)
        if (!(mapped_iti %in% trial_names)) {
          stop(
            paste(
              "Contingency", i, "references unknown ITI trial:",
              raw_parts[8]
            )
          )
        }
        raw_parts[8] <- mapped_iti
      }

      paste(raw_parts, collapse = ", ")
    }, character(1))
  }

  get_phase_order <- function(phase_values) {
    phase_values <- as.character(phase_values)
    if (length(phase_values) == 0) {
      return(character(0))
    }

    cont_cfg <- contingencies()
    if ((is.null(cont_cfg) || length(cont_cfg) == 0) && exists("contingencies", envir = .GlobalEnv)) {
      cont_cfg <- get("contingencies", envir = .GlobalEnv)
    }

    if (!is.null(cont_cfg) && length(cont_cfg) > 0) {
      configured_phases <- unique(vapply(cont_cfg, function(x) {
        trimws(unlist(strsplit(x, ","))[1])
      }, character(1)))

      ordered <- configured_phases[configured_phases %in% phase_values]
      remaining <- setdiff(unique(phase_values), ordered)
      return(c(ordered, remaining))
    }

    unique(phase_values)
  }

  workflow_tabs <- c("home", "network", "trials", "contingencies", "simulate", "results_individual", "results_general", "help")
  workflow_labels <- c(
    home = "Home",
    network = "Network",
    trials = "Trials",
    contingencies = "Contingencies",
    simulate = "Simulate",
    results_individual = "Individual",
    results_general = "General",
    help = "Help"
  )

  screen_hints <- list(
    home = list(
      title = "Start Here",
      text = "Use the minimal extinction template for a guided first run, then move to Network Architecture to inspect or edit."
    ),
    network = list(
      title = "Build the Architecture",
      text = "Define NPUs and connections first. Use the network canvas to validate structure before pressing Create Architecture."
    ),
    trials = list(
      title = "Define Trial Timesteps",
      text = "Set the number of timesteps, define sensory/US values per timestep, and save trial types for contingencies."
    ),
    contingencies = list(
      title = "Configure Phases",
      text = "Add phases with trial presentation rules. Enable ITI only when needed for your experimental design."
    ),
    simulate = list(
      title = "Run and Persist",
      text = "Run one or multiple networks, then save complete sessions to RDS to preserve architecture, contingencies, and outputs."
    ),
    results_individual = list(
      title = "Inspect a Single Network",
      text = "Plot activations and weights for one network to diagnose learning trajectories and extinction shape."
    ),
    results_general = list(
      title = "Compare Across Networks",
      text = "Use aggregated plots to check variability and central tendency across all simulated networks."
    ),
    help = list(
      title = "Reference Guide",
      text = "Use this section as a glossary of fields, file formats, and recommended workflow patterns."
    )
  )

  output$workflow_progress <- renderUI({
    current_tab <- input$main_nav
    if (is.null(current_tab) || !current_tab %in% workflow_tabs) {
      current_tab <- "home"
    }
    active_idx <- match(current_tab, workflow_tabs)

    tags$div(
      class = "ux-flow",
      lapply(seq_along(workflow_tabs), function(i) {
        state_class <- if (i < active_idx) {
          "done"
        } else if (i == active_idx) {
          "active"
        } else {
          ""
        }
        tags$div(
          class = trimws(paste("ux-step", state_class)),
          tags$span(class = "ux-step-num", i),
          tags$span(class = "ux-step-label", workflow_labels[[workflow_tabs[i]]])
        )
      })
    )
  })

  output$screen_hint <- renderUI({
    current_tab <- input$main_nav
    if (is.null(current_tab) || !current_tab %in% names(screen_hints)) {
      current_tab <- "home"
    }
    hint <- screen_hints[[current_tab]]
    mode_note <- if (current_tab == "network" && app_mode() == "Beginner") {
      "Beginner mode keeps appendix defaults for threshold and learning-rate parameters."
    } else {
      ""
    }
    hint_text <- paste(c(hint$text, mode_note), collapse = " ")
    hint_text <- gsub("\\s+", " ", trimws(hint_text))

    tags$div(
      class = "ux-hint",
      tags$span(class = "ux-hint-icon", tags$i(class = "fa fa-lightbulb-o")),
      tags$div(
        tags$div(class = "ux-hint-title", hint$title),
        tags$p(class = "ux-hint-text", hint_text)
      )
    )
  })

  # Función para ordenar unidades
  get_unit_order <- function(units, npes_data) {
    layer_order <- c("PrimarySensory", "AssociativeSensory", "AssociativeMotor", "PrimaryMotor", "Hippocampal", "Dopaminergic", "US")

    # Crear un data frame con unidades y sus capas correspondientes
    unit_layers <- npes_data %>%
      select(NPE, Layer) %>%
      filter(NPE %in% units)

    # Ordenar las unidades según el orden de capas definido
    ordered_units <- unit_layers %>%
      mutate(LayerOrder = match(Layer, layer_order)) %>%
      arrange(LayerOrder, NPE) %>%
      pull(NPE)

    # Añadir cualquier unidad que no esté en npes_data al final
    remaining_units <- setdiff(units, ordered_units)
    c(ordered_units, sort(remaining_units))
  }

  auto_layout_by_layer <- function(npe_df) {
    if (nrow(npe_df) == 0) {
      return(npe_df)
    }

    if (!"x" %in% names(npe_df)) npe_df$x <- NA_real_
    if (!"y" %in% names(npe_df)) npe_df$y <- NA_real_

    layer_x <- c(
      PrimarySensory = -80,
      AssociativeSensory = 180,
      Hippocampal = 300,
      AssociativeMotor = 460,
      PrimaryMotor = 650,
      Dopaminergic = 460,
      US = -260
    )

    layer_center_y <- c(
      PrimarySensory = -20,
      AssociativeSensory = -20,
      Hippocampal = 160,
      AssociativeMotor = -20,
      PrimaryMotor = -20,
      Dopaminergic = 250,
      US = 330
    )

    laid_out <- npe_df %>%
      arrange(Layer, NPE) %>%
      group_by(Layer) %>%
      mutate(
        idx = row_number(),
        n_layer = n(),
        y = unname(layer_center_y[Layer]) + (idx - (n_layer + 1) / 2) * 105,
        x = unname(layer_x[Layer])
      ) %>%
      ungroup()

    laid_out %>%
      select(-idx, -n_layer)
  }

  get_extinction_template <- function() {
    list(
      NPEs = data.frame(
        NPE = c("US", "D", "S1", "S..1", "H1", "M..1", "M.1"),
        Type = rep("Excitatory", 7),
        Layer = c("US", "Dopaminergic", "PrimarySensory", "AssociativeSensory", "Hippocampal", "AssociativeMotor", "PrimaryMotor"),
        Activation = rep(0, 7),
        Temporal.Summation = rep(0.1, 7),
        Activation.Decay = rep(0.1, 7),
        mu = rep(0.2, 7),
        sigma = rep(0.15, 7),
        logisSigma = rep(0.1, 7),
        x = rep(NA_real_, 7),
        y = rep(NA_real_, 7),
        stringsAsFactors = FALSE
      ),
      Connections = data.frame(
        PreSinapticNPE = c("S1", "S..1", "S..1", "M..1", "M..1", "US"),
        PostSinapticNPE = c("S..1", "H1", "M..1", "D", "M.1", "D"),
        Weight = c(0.1, 0.1, 0.1, 0.1, 0.1, 1.0),
        alpha = rep(0.5, 6),
        beta = rep(0.12, 6),
        alpha_prime = rep(0.5, 6),
        beta_prime = rep(0.12, 6),
        stringsAsFactors = FALSE
      ),
      trials = list(
        Training = c(
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True",
          "US,1.00,S1,1.00,True"
        ),
        Extinction = c(
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True"
        )
      ),
      contingencies = c(
        "training, Random, Training, 100, False",
        "extinction, Random, Extinction, 100, False"
      )
    )
  }

  observeEvent(input$load_extinction_template, {
    template <- get_extinction_template()
    template_npes <- auto_layout_by_layer(template$NPEs)
    npes(template_npes)
    connections(template$Connections)
    trials(template$trials)
    contingencies(template$contingencies)
    simulation_results(NULL)

    assign("NPEs", template_npes, envir = .GlobalEnv)
    assign("Connections", template$Connections, envir = .GlobalEnv)
    assign("trials", template$trials, envir = .GlobalEnv)
    assign("contingencies", template$contingencies, envir = .GlobalEnv)
    if (exists("datos_simulacion", envir = .GlobalEnv)) rm("datos_simulacion", envir = .GlobalEnv)

    updateSelectInput(session, "conn_pre", choices = template_npes$NPE)
    updateSelectInput(session, "conn_post", choices = template_npes$NPE)
    updateSelectInput(session, "trial_types", choices = names(template$trials))
    updateSelectInput(session, "iti_trial", choices = character(0))
    updateSelectInput(session, "selected_simulation", choices = character(0))

    showNotification("Minimal extinction template loaded. You can run the simulation directly.", type = "message")
  })

  observeEvent(TRUE, {
    npes(auto_layout_by_layer(npes()))
    set_all_path_inputs(current_workdir())
    apply_mode(app_mode())
  }, once = TRUE)

  observeEvent(input$app_mode, {
    app_mode(input$app_mode)
    apply_mode(input$app_mode)
    showNotification(paste("Mode changed to", input$app_mode), type = "message")
  }, ignoreInit = TRUE)

  observeEvent(input$use_workdir_paths, {
    set_all_path_inputs(current_workdir())
    showNotification("All path fields now use the current folder.", type = "message")
  })

  observeEvent(input$use_downloads_paths, {
    set_all_path_inputs(downloads_dir())
    showNotification("All path fields now use Downloads.", type = "message")
  })

  # Añadir NPE
  observeEvent(input$add_npe, {
    npe_mu <- if (app_mode() == "Beginner") 0.2 else input$npe_mu
    npe_sigma <- if (app_mode() == "Beginner") 0.15 else input$npe_sigma
    npe_logis <- if (app_mode() == "Beginner") 0.1 else input$npe_logistic_slope

    new_npe <- data.frame(
      NPE = input$npe_name,
      Type = input$npe_type,
      Layer = input$npe_layer,
      Activation = input$npe_activation,
      Temporal.Summation = input$npe_temporal_summation,
      Activation.Decay = input$npe_decay,
      mu = npe_mu,
      sigma = npe_sigma,
      logisSigma = npe_logis,
      x = NA_real_,
      y = NA_real_,
      stringsAsFactors = FALSE
    )
    current_npes <- auto_layout_by_layer(rbind(npes(), new_npe))
    npes(current_npes)
    updateSelectInput(session, "conn_pre", choices = current_npes$NPE)
    updateSelectInput(session, "conn_post", choices = current_npes$NPE)
  })

  # Añadir Conexión
  observeEvent(input$add_connection, {
    req(input$conn_pre != input$conn_post)

    conn_alpha <- if (app_mode() == "Beginner") 0.5 else input$conn_alpha
    conn_beta <- if (app_mode() == "Beginner") 0.12 else input$conn_beta
    conn_alpha_prime <- if (app_mode() == "Beginner") 0.5 else input$conn_alpha_prime
    conn_beta_prime <- if (app_mode() == "Beginner") 0.12 else input$conn_beta_prime

    new_connection <- data.frame(
      PreSinapticNPE = input$conn_pre,
      PostSinapticNPE = input$conn_post,
      Weight = input$conn_weight,
      alpha = conn_alpha,
      beta = conn_beta,
      alpha_prime = conn_alpha_prime,
      beta_prime = conn_beta_prime,
      stringsAsFactors = FALSE
    )
    current_connections <- rbind(connections(), new_connection)
    connections(current_connections)
  })

  # Tabla de NPEs
  output$npe_table <- renderDT({
    datatable(npes(), options = list(pageLength = 5, scrollX = TRUE, scrollY = "200px"))
  })

  # Tabla de Conexiones
  output$connection_table <- renderDT({
    datatable(connections(), options = list(pageLength = 5, scrollX = TRUE, scrollY = "200px"))
  })

  # Visualización de la red
  output$network_plot <- renderVisNetwork({
    req(nrow(npes()) > 0, nrow(connections()) > 0)

    base_nodes <- npes()
    if (
      !all(c("x", "y") %in% names(base_nodes)) ||
      any(!is.finite(base_nodes$x)) ||
      any(!is.finite(base_nodes$y))
    ) {
      base_nodes <- auto_layout_by_layer(base_nodes)
      npes(base_nodes)
    }

    inhibitory_units <- base_nodes %>%
      filter(Type == "Inhibitory") %>%
      pull(NPE)

    nodes <- base_nodes %>%
      mutate(
        id = NPE,
        label = NPE,
        shape = case_when(
          Layer == "US" ~ "hexagon",
          Layer == "PrimarySensory" ~ "square",
          Type == "Inhibitory" ~ "diamond",
          TRUE ~ "circle"
        ),
        color.background = "#ffffff",
        color.border = case_when(
          Layer %in% c("PrimarySensory", "US") ~ "#111111",
          Layer == "Dopaminergic" ~ "#111111",
          Type == "Inhibitory" ~ "#111111",
          TRUE ~ "#111111"
        ),
        borderWidth = case_when(
          Layer == "US" ~ 3.5,
          Layer == "Dopaminergic" ~ 3,
          Type == "Inhibitory" ~ 2.5,
          TRUE ~ 2
        ),
        size = case_when(
          Layer %in% c("PrimarySensory", "US") ~ 22,
          Type == "Inhibitory" ~ 26,
          TRUE ~ 30
        ),
        x = as.numeric(x),
        y = as.numeric(y),
        fixed = FALSE,
        title = paste(
          "<b>NPU:</b>", NPE,
          "<br><b>Type:</b>", Type,
          "<br><b>Layer:</b>", Layer
        )
      )

    edges <- connections() %>%
      mutate(
        from = PreSinapticNPE,
        to = PostSinapticNPE,
        inhibitory = PreSinapticNPE %in% inhibitory_units,
        fixed_link = Weight >= 0.999,
        arrows = "to",
        color = case_when(
          fixed_link ~ "#111111",
          inhibitory ~ "#111111",
          TRUE ~ "#666666"
        ),
        dashes = !fixed_link & !inhibitory,
        width = ifelse(fixed_link, 5, ifelse(inhibitory, 2.5, pmax(1.5, Weight * 3.5))),
        smooth = ifelse(inhibitory, "curvedCW", "continuous"),
        title = paste(
          "<b>Connection:</b>", PreSinapticNPE, "→", PostSinapticNPE,
          "<br><b>Weight:</b>", round(Weight, 4),
          "<br><b>alpha:</b>", alpha,
          "<br><b>beta:</b>", beta
        )
      )

    visNetwork(nodes, edges, width = "100%", height = "520px") %>%
      visNodes(
        shadow = list(enabled = FALSE),
        font = list(size = 28, face = "bold", color = "#111111")
      ) %>%
      visEdges(arrows = list(to = list(enabled = TRUE, scaleFactor = 0.75))) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = TRUE) %>%
      visPhysics(enabled = FALSE) %>%
      visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE, navigationButtons = FALSE, hover = TRUE) %>%
      visExport(type = "png", name = "ddm_architecture", background = "#ffffff") %>%
      visEvents(type = "once", startStabilizing = "function() {
          this.moveTo({scale:1.0})
        }")
  })

  # Actualizar posiciones de los nodos
  observe({
    network_data <- input$network_plot_positions
    if (!is.null(network_data) && !is.null(network_data$nodes)) {
      node_positions <- as.data.frame(network_data$nodes)
      required_cols <- c("id", "x", "y")
      if (!all(required_cols %in% names(node_positions))) {
        return()
      }
      updated_npes <- npes()
      idx <- match(updated_npes$NPE, node_positions$id)
      valid <- !is.na(idx)
      updated_npes$x[valid] <- as.numeric(node_positions$x[idx[valid]])
      updated_npes$y[valid] <- as.numeric(node_positions$y[idx[valid]])
      npes(updated_npes)
    }
  })

  # Reorganizar la red
  observeEvent(input$reorganize_network, {
    npes(auto_layout_by_layer(npes()))
  })

  # Crear Arquitectura
  observeEvent(input$create_architecture, {
    req(nrow(npes()) > 0, nrow(connections()) > 0)

    # Guardar los data frames y la lista de ensayos
    assign("NPEs", isolate(npes()), envir = .GlobalEnv)
    assign("Connections", isolate(connections()), envir = .GlobalEnv)
    assign("trials", isolate(trials()), envir = .GlobalEnv)

    showNotification("Architecture created successfully", type = "message")
  })

  # Guardar NPEs
  observeEvent(input$save_npes, {
    req(input$npe_file_name != "")
    tryCatch(
      {
        filename <- ensure_csv_extension(input$npe_file_name)
        full_path <- file.path(resolve_directory(input$npe_file_path), filename)
        write.csv(npes(), full_path, row.names = FALSE)
        showNotification("NPEs saved successfully", type = "message")
      },
      error = function(e) {
        showNotification(paste("Error saving NPEs:", e$message), type = "error", duration = NULL)
      }
    )
  })

  # Guardar Conexiones
  observeEvent(input$save_connections, {
    req(input$conn_file_name != "")
    tryCatch(
      {
        filename <- ensure_csv_extension(input$conn_file_name)
        full_path <- file.path(resolve_directory(input$conn_file_path), filename)
        write.csv(connections(), full_path, row.names = FALSE)
        showNotification("Connections saved successfully", type = "message")
      },
      error = function(e) {
        showNotification(paste("Error saving Connections:", e$message), type = "error", duration = NULL)
      }
    )
  })

  # Importar NPEs
  observeEvent(input$import_npes, {
    req(input$npe_file_name != "")
    tryCatch(
      {
        filename <- ensure_csv_extension(input$npe_file_name)
        full_path <- file.path(resolve_directory(input$npe_file_path), filename)
        if (!file.exists(full_path)) {
          stop(paste("The file does not exist:", full_path))
        }
        imported_npes <- as.data.frame(read_csv(full_path))
        if (!"x" %in% names(imported_npes)) imported_npes$x <- NA_real_
        if (!"y" %in% names(imported_npes)) imported_npes$y <- NA_real_
        imported_npes <- auto_layout_by_layer(imported_npes)
        npes(imported_npes)
        updateSelectInput(session, "conn_pre", choices = imported_npes$NPE)
        updateSelectInput(session, "conn_post", choices = imported_npes$NPE)
        showNotification("NPEs imported successfully", type = "message")
      },
      error = function(e) {
        showNotification(paste("Error importing NPEs:", e$message), type = "error", duration = NULL)
      }
    )
  })

  # Importar Conexiones
  observeEvent(input$import_connections, {
    req(input$conn_file_name != "")
    tryCatch(
      {
        filename <- ensure_csv_extension(input$conn_file_name)
        full_path <- file.path(resolve_directory(input$conn_file_path), filename)
        if (!file.exists(full_path)) {
          stop(paste("The file does not exist:", full_path))
        }
        imported_connections <- as.data.frame(read_csv(full_path))
        connections(imported_connections)
        showNotification("Connections imported successfully", type = "message")
      },
      error = function(e) {
        showNotification(paste("Error importing Connections:", e$message), type = "error", duration = NULL)
      }
    )
  })

  # Observador para el número de momentos temporales
  observeEvent(input$num_moments,
    {
      if (is.null(input$num_moments) || is.na(input$num_moments) || input$num_moments < 1) {
        updateNumericInput(session, "num_moments", value = 1)
        showNotification("The minimum number of time steps is 1", type = "warning")
      }
    },
    ignoreInit = TRUE,
    ignoreNULL = FALSE
  )

  # Observadores para limitar los valores de entrada a 0-1
  observe({
    req(input$num_moments)
    lapply(1:input$num_moments, function(moment) {
      lapply(npes()$NPE[npes()$Layer %in% c("US", "PrimarySensory")], function(npe) {
        input_id <- paste0("input_", npe, "_", moment)
        observeEvent(input[[input_id]],
          {
            if (is.null(input[[input_id]]) || is.na(input[[input_id]]) || input[[input_id]] < 0) {
              updateNumericInput(session, input_id, value = 0)
              showNotification("Values must be between 0 and 1", type = "warning")
            } else if (input[[input_id]] > 1) {
              updateNumericInput(session, input_id, value = 1)
              showNotification("Values must be between 0 and 1", type = "warning")
            }
          },
          ignoreInit = TRUE,
          ignoreNULL = FALSE
        )
      })
    })
  })

  # Regenerar inputs dinámicos cuando cambia el número de momentos
  observeEvent(input$num_moments,
    {
      output$dynamic_inputs <- renderUI({
        req(npes(), input$num_moments)
        input_list <- list()
        valid_npes <- npes()$NPE[npes()$Layer %in% c("US", "PrimarySensory")]

        for (moment in 1:input$num_moments) {
          input_list[[paste0("moment_", moment)]] <- list(
            h4(paste("Timestep", moment)),
            fluidRow(
              lapply(valid_npes, function(npe) {
                column(
                  width = 12 / length(valid_npes),
                  numericInput(
                    paste0("input_", npe, "_", moment),
                    label = npe,
                    value = 0,
                    min = 0,
                    max = 1,
                    step = 0.1
                  )
                )
              })
            )
          )
        }

        input_list
      })
    },
    ignoreInit = TRUE
  )

  # Función para crear la cadena del ensayo
  create_trial_string <- function(input_values, learning_rule) {
    trial_string <- paste(names(input_values), sprintf("%.2f", input_values), sep = ",", collapse = ",")
    paste0(trial_string, ",", ifelse(learning_rule, "True", "False"))
  }

  # Añadir ensayo
  observeEvent(input$add_trial, {
    req(input$trial_name)
    valid_npes <- npes()$NPE[npes()$Layer %in% c("US", "PrimarySensory")]
    new_trials <- sapply(1:input$num_moments, function(moment) {
      input_values <- sapply(valid_npes, function(npe) {
        value <- input[[paste0("input_", npe, "_", moment)]]
        if (is.null(value) || is.na(value)) 0 else min(value, 1)
      })
      create_trial_string(input_values, input$learning_rule)
    })
    current_trials <- isolate(trials())
    current_trials[[input$trial_name]] <- new_trials
    trials(current_trials)
    assign("trials", isolate(trials()), envir = .GlobalEnv)
    showNotification("Trial added successfully", type = "message")

    # Actualizar las opciones de tipos de ensayos en la pestaña de contingencias
    updateSelectInput(session, "trial_types", choices = names(current_trials))
  })

  # Añadir ITI
  observeEvent(input$add_iti, {
    showModal(modalDialog(
      title = "Configure ITI",
      textInput("iti_name", "ITI Name"),
      uiOutput("iti_inputs"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_iti", "Save ITI", icon = icon("save"))
      )
    ))
  })

  # Generar inputs dinámicos para el ITI
  output$iti_inputs <- renderUI({
    req(npes())
    input_list <- list()
    valid_npes <- npes()$NPE[npes()$Layer %in% c("US", "PrimarySensory")]
    for (npe in valid_npes) {
      input_list[[npe]] <- numericInput(paste0("iti_input_", npe),
        label = npe,
        value = 0,
        min = 0,
        max = 1,
        step = 0.1
      )
    }
    input_list$iti_learning_rule <- checkboxInput("iti_learning_rule", "Active learning rule", value = TRUE)
    input_list
  })

  # Observador para los inputs del ITI
  observe({
    req(input$iti_name)
    valid_npes <- npes()$NPE[npes()$Layer %in% c("US", "PrimarySensory")]
    lapply(valid_npes, function(npe) {
      input_id <- paste0("iti_input_", npe)
      observeEvent(input[[input_id]], {
        if (is.na(input[[input_id]]) || input[[input_id]] < 0) {
          updateNumericInput(session, input_id, value = 0)
          showNotification("Values must be between 0 and 1", type = "warning")
        } else if (input[[input_id]] > 1) {
          updateNumericInput(session, input_id, value = 1)
          showNotification("Values must be between 0 and 1", type = "warning")
        }
      })
    })
  })

  # Guardar ITI
  observeEvent(input$save_iti, {
    req(input$iti_name)
    valid_npes <- npes()$NPE[npes()$Layer %in% c("US", "PrimarySensory")]
    input_values <- sapply(valid_npes, function(npe) {
      value <- input[[paste0("iti_input_", npe)]]
      if (is.null(value) || is.na(value)) 0 else min(value, 1)
    })
    iti_trial <- create_trial_string(input_values, input$iti_learning_rule)
    current_trials <- isolate(trials())
    current_trials[[input$iti_name]] <- iti_trial
    trials(current_trials)
    assign("trials", isolate(trials()), envir = .GlobalEnv)
    removeModal()
    showNotification("ITI added successfully", type = "message")

    # Actualizar las opciones de tipos de ensayos en la pestaña de contingencias
    updateSelectInput(session, "trial_types", choices = names(current_trials))
    updateSelectInput(session, "iti_trial",
      choices = names(current_trials)[sapply(current_trials, function(x) is.character(x) && length(x) == 1)]
    )
  })

  # Mostrar botones de ensayos
  output$trial_buttons <- renderUI({
    all_trials <- trials()
    if (length(all_trials) == 0) {
      return(NULL)
    }

    tagList(
      lapply(names(all_trials), function(trial_type) {
        div(
          actionButton(
            inputId = paste0("trial_", trial_type),
            label = trial_type,
            class = "btn-block"
          ),
          div(id = paste0("content_", trial_type), style = "display: none;")
        )
      })
    )
  })

  # Función auxiliar para renderizar el contenido del ensayo
  renderTrialContent <- function(trial_type, trial_data) {
    if (is.character(trial_data) && length(trial_data) == 1) {
      # Para ITI
      trial_parts <- unlist(strsplit(trial_data, ","))
      content <- tags$div(
        h5(paste("ITI details:", trial_type)),
        lapply(seq(1, length(trial_parts) - 1, 2), function(i) {
          p(paste(trial_parts[i], ":", trial_parts[i + 1]))
        }),
        p(paste("Learning rule:", trial_parts[length(trial_parts)]))
      )
    } else {
      # Para otros tipos de ensayos
      content <- tags$div(
        h5(paste("Trial details:", trial_type)),
        lapply(seq_along(trial_data), function(i) {
          trial_parts <- unlist(strsplit(trial_data[i], ","))
          tags$div(
            h6(paste("Timestep", i)),
            lapply(seq(1, length(trial_parts) - 1, 2), function(j) {
              p(paste(trial_parts[j], ":", trial_parts[j + 1]))
            }),
            p(paste("Learning rule:", trial_parts[length(trial_parts)]))
          )
        })
      )
    }
    as.character(content)
  }

  # Manejar clics en los botones de ensayos
  observe({
    all_trials <- trials()
    lapply(names(all_trials), function(trial_type) {
      observeEvent(input[[paste0("trial_", trial_type)]], {
        content_id <- paste0("content_", trial_type)
        if (is.null(input[[paste0("trial_", trial_type)]]) || input[[paste0("trial_", trial_type)]] %% 2 == 1) {
          # Mostrar contenido
          shinyjs::show(content_id)
          # Actualizar el contenido
          trial_data <- all_trials[[trial_type]]
          content <- renderTrialContent(trial_type, trial_data)
          shinyjs::html(content_id, content)
        } else {
          # Ocultar contenido
          shinyjs::hide(content_id)
        }
      })
    })
  })

  # Generar inputs dinámicos para el número de ensayos
  output$trial_numbers <- renderUI({
    req(input$trial_types)
    lapply(input$trial_types, function(trial_type) {
      numericInput(paste0("num_", trial_type),
        label = paste("Number of trials from", trial_type),
        value = 1, min = 1
      )
    })
  })

  # Añadir Contingencia
  observeEvent(input$add_contingency, {
    # Verificar que todos los campos necesarios estén llenos
    if (input$phase_name == "") {
      showNotification("The phase or condition name cannot be left blank", type = "error")
      return()
    }

    if (length(input$trial_types) == 0) {
      showNotification("You must select at least one trial type", type = "error")
      return()
    }

    # Verificar que se haya ingresado un número de ensayos para cada tipo de ensayo seleccionado
    for (tt in input$trial_types) {
      if (is.null(input[[paste0("num_", tt)]]) || is.na(input[[paste0("num_", tt)]]) || input[[paste0("num_", tt)]] < 1) {
        showNotification(paste("The number of trials for", tt, "cannot be left blank or less than 1"), type = "error")
        return()
      }
    }

    # Si no se resetean las activaciones, verificar los campos adicionales
    if (!input$reset_activations) {
      if (is.null(input$min_iti) || is.na(input$min_iti) || input$min_iti < 1) {
        showNotification("The minimum ITI value cannot be left blank or less than 1", type = "error")
        return()
      }
      if (is.null(input$max_iti) || is.na(input$max_iti) || input$max_iti < 1) {
        showNotification("The maximum ITI value cannot be left blank or less than 1", type = "error")
        return()
      }
      if (is.null(input$iti_trial) || input$iti_trial == "") {
        showNotification("You must select a created ITI", type = "error")
        return()
      }
    }

    # Si todas las validaciones pasan, proceder con la creación de la contingencia
    trial_types <- paste(input$trial_types, collapse = "/")
    trial_numbers <- paste(sapply(input$trial_types, function(tt) input[[paste0("num_", tt)]]), collapse = "-")

    reset_activations <- input$reset_activations

    if (!reset_activations) {
      contingency <- paste(
        input$phase_name,
        input$presentation_mode,
        trial_types,
        trial_numbers,
        "True",
        input$min_iti,
        input$max_iti,
        input$iti_trial,
        sep = ", "
      )
    } else {
      contingency <- paste(
        input$phase_name,
        input$presentation_mode,
        trial_types,
        trial_numbers,
        "False",
        sep = ", "
      )
    }

    current_contingencies <- c(isolate(contingencies()), contingency)
    contingencies(current_contingencies)
    assign("contingencies", current_contingencies, envir = .GlobalEnv)

    showNotification("Contingency added successfully", type = "message")
  })

  # Mostrar tabla de contingencies
  output$contingencies_table <- renderDT({
    cont <- contingencies()
    if (length(cont) == 0) {
      return(NULL)
    }

    cont_df <- data.frame(
      Contingency = cont,
      stringsAsFactors = FALSE
    )

    datatable(cont_df, options = list(pageLength = 10, scrollX = TRUE, scrollY = "300px"))
  })

  # Observador para el número de simulaciones
  observeEvent(input$num_simulations,
    {
      if (is.na(input$num_simulations) || input$num_simulations < 1) {
        updateNumericInput(session, "num_simulations", value = 1)
        showNotification("The minimum number of simulations is 1", type = "warning")
      }
    },
    ignoreInit = TRUE
  )

  # Ejecutar simulación
  observeEvent(input$run_simulation, {
    # Verificar que todos los datos necesarios existan
    missing_components <- character(0)

    if (!exists("NPEs", envir = .GlobalEnv) || nrow(get("NPEs", envir = .GlobalEnv)) == 0) {
      missing_components <- c(missing_components, "NPEs")
    }
    if (!exists("Connections", envir = .GlobalEnv) || nrow(get("Connections", envir = .GlobalEnv)) == 0) {
      missing_components <- c(missing_components, "Connections")
    }
    if (!exists("trials", envir = .GlobalEnv) || length(get("trials", envir = .GlobalEnv)) == 0) {
      missing_components <- c(missing_components, "trials")
    }
    if (!exists("contingencies", envir = .GlobalEnv) || length(get("contingencies", envir = .GlobalEnv)) == 0) {
      missing_components <- c(missing_components, "contingencies")
    }

    if (length(missing_components) > 0) {
      showNotification(
        paste(
          "Missing necessary data for the simulation:",
          paste(missing_components, collapse = ", "),
          ". Ensure that you have created the architecture and defined the trials and contingencies."
        ),
        type = "error", duration = NULL
      )
      return()
    }

    output$simulation_status <- renderText("Simulation in progress...")
    updateProgressBar(session, "sim_progress", value = 0)

    # Preparación previa compartida por todas las redes
    NPEs <- as.data.frame(get("NPEs", envir = .GlobalEnv), stringsAsFactors = FALSE)
    Connections <- as.data.frame(get("Connections", envir = .GlobalEnv), stringsAsFactors = FALSE)
    trial_list <- get("trials", envir = .GlobalEnv)
    contingency_list <- get("contingencies", envir = .GlobalEnv)

    missing_npes <- verify_npes(trial_list, NPEs)
    if (length(missing_npes) > 0) {
      missing_report <- paste(
        vapply(names(missing_npes), function(tt) {
          paste0(tt, ": ", paste(unique(missing_npes[[tt]]), collapse = "/"))
        }, character(1)),
        collapse = "; "
      )
      msg <- paste("Trials include NPUs not present in architecture:", missing_report)
      output$simulation_status <- renderText(msg)
      showNotification(msg, type = "error", duration = NULL)
      return()
    }

    prep <- tryCatch(
      {
        normalized_contingencies <- harmonize_contingencies(contingency_list, trial_list)
        HasITI <- sapply(normalized_contingencies, function(x) {
          as.logical(trimws(unlist(strsplit(x, ",")))[5])
        })
        TimeSteps <- Create.Phases(normalized_contingencies, trial_list)
        list(
          contingencies = normalized_contingencies,
          HasITI = HasITI,
          TimeSteps = TimeSteps
        )
      },
      error = function(e) {
        msg <- paste("Error while preparing contingencies:", e$message)
        output$simulation_status <- renderText(msg)
        showNotification(msg, type = "error", duration = NULL)
        NULL
      }
    )

    if (is.null(prep)) {
      return()
    }

    assign("contingencies", prep$contingencies, envir = .GlobalEnv)
    assign("HasITI", prep$HasITI, envir = .GlobalEnv)
    assign("TimeSteps", prep$TimeSteps, envir = .GlobalEnv)

    all_results <- vector("list", input$num_simulations)
    failed_networks <- integer(0)
    failed_messages <- character(0)

    for (i in 1:input$num_simulations) {
      tryCatch(
        {
          datos <- Simulate.DBP(NPEs, Connections, prep$TimeSteps, prep$HasITI)
          all_results[[i]] <- datos

          updateProgressBar(session, "sim_progress", value = (i / input$num_simulations) * 100)
          output$simulation_status <- renderText(paste("Simulation", i, "of", input$num_simulations, "completed"))
        },
        error = function(e) {
          failed_networks <<- c(failed_networks, i)
          failed_messages <<- c(failed_messages, paste0("Network ", i, ": ", e$message))
          output$simulation_status <- renderText(paste("Error in simulation", i, ":", e$message))
        }
      )
    }

    successful_results <- Filter(function(x) !is.null(x), all_results)
    if (length(successful_results) == 0) {
      simulation_results(NULL)
      if (exists("datos_simulacion", envir = .GlobalEnv)) {
        rm("datos_simulacion", envir = .GlobalEnv)
      }
      updateSelectInput(session, "selected_simulation", choices = character(0), selected = character(0))

      detail <- if (length(failed_messages) > 0) failed_messages[[1]] else "Unknown simulation error."
      output$simulation_status <- renderText(paste("Simulation failed for all networks.", detail))
      showNotification(
        paste("No network finished successfully.", detail),
        type = "error",
        duration = NULL
      )
      return()
    }

    # Guardar todos los resultados válidos
    simulation_results(successful_results)
    assign("datos_simulacion", successful_results, envir = .GlobalEnv)

    network_choices <- paste("Network", seq_along(successful_results))
    updateSelectInput(session, "selected_simulation",
      choices = network_choices,
      selected = network_choices[1]
    )

    if (length(failed_networks) > 0) {
      output$simulation_status <- renderText(
        paste0("Completed ", length(successful_results), " of ", input$num_simulations, " networks")
      )
      showNotification(
        paste("Networks with errors:", paste(failed_networks, collapse = ",")),
        type = "warning",
        duration = NULL
      )
    } else {
      output$simulation_status <- renderText("All simulations completed")
    }
  })

  # Guardar simulación
  observeEvent(input$save_simulation, {
    req(input$sim_file_name != "")

    missing_components <- character(0)

    if (!exists("NPEs", envir = .GlobalEnv) || nrow(get("NPEs", envir = .GlobalEnv)) == 0) {
      missing_components <- c(missing_components, "NPEs")
    }
    if (!exists("Connections", envir = .GlobalEnv) || nrow(get("Connections", envir = .GlobalEnv)) == 0) {
      missing_components <- c(missing_components, "Connections")
    }
    if (!exists("trials", envir = .GlobalEnv) || length(get("trials", envir = .GlobalEnv)) == 0) {
      missing_components <- c(missing_components, "trials")
    }
    if (!exists("contingencies", envir = .GlobalEnv) || length(get("contingencies", envir = .GlobalEnv)) == 0) {
      missing_components <- c(missing_components, "contingencies")
    }
    if (!exists("datos_simulacion", envir = .GlobalEnv) || length(get("datos_simulacion", envir = .GlobalEnv)) == 0) {
      missing_components <- c(missing_components, "Simulation data")
    }

    if (length(missing_components) > 0) {
      showNotification(
        paste(
          "The simulation cannot be saved. The following components are missing:",
          paste(missing_components, collapse = ", ")
        ),
        type = "error", duration = NULL
      )
      return()
    }

    simulation_data <- list(
      NPEs = get("NPEs", envir = .GlobalEnv),
      Connections = get("Connections", envir = .GlobalEnv),
      trials = get("trials", envir = .GlobalEnv),
      contingencies = get("contingencies", envir = .GlobalEnv),
      datos_simulacion = get("datos_simulacion", envir = .GlobalEnv)
    )

    tryCatch(
      {
        filename <- ensure_rds_extension(input$sim_file_name)
        full_path <- file.path(resolve_directory(input$sim_file_path), filename)
        saveRDS(simulation_data, file = full_path)
        showNotification("Simulation saved successfully", type = "message")
      },
      error = function(e) {
        showNotification(paste("Error saving the simulation:", e$message), type = "error", duration = NULL)
      }
    )
  })

  # Cargar simulación
  observeEvent(input$load_simulation, {
    req(input$sim_file_name != "")
    tryCatch(
      {
        filename <- ensure_rds_extension(input$sim_file_name)
        full_path <- file.path(resolve_directory(input$sim_file_path), filename)
        if (!file.exists(full_path)) {
          stop(paste("The file does not exist:", full_path))
        }
        simulation_data <- readRDS(full_path)

        # Actualizar los valores reactivos y el entorno global
        npes(simulation_data$NPEs)
        connections(simulation_data$Connections)
        trials(simulation_data$trials)
        contingencies(simulation_data$contingencies)
        simulation_results(simulation_data$datos_simulacion)

        assign("NPEs", simulation_data$NPEs, envir = .GlobalEnv)
        assign("Connections", simulation_data$Connections, envir = .GlobalEnv)
        assign("trials", simulation_data$trials, envir = .GlobalEnv)
        assign("contingencies", simulation_data$contingencies, envir = .GlobalEnv)
        assign("datos_simulacion", simulation_data$datos_simulacion, envir = .GlobalEnv)

        # Actualizar las opciones de tipos de ensayos en la pestaña de contingencias
        updateSelectInput(session, "trial_types", choices = names(simulation_data$trials))
        updateSelectInput(session, "iti_trial",
          choices = names(simulation_data$trials)[sapply(simulation_data$trials, function(x) is.character(x) && length(x) == 1)]
        )

        # Actualizar el selector de simulaciones
        loaded_choices <- paste("Network", seq_along(simulation_data$datos_simulacion))
        updateSelectInput(session, "selected_simulation",
          choices = loaded_choices,
          selected = if (length(loaded_choices) > 0) loaded_choices[1] else character(0)
        )

        showNotification("Simulation successfully loaded", type = "message")
      },
      error = function(e) {
        showNotification(paste("Error loading the simulation:", e$message), type = "error", duration = NULL)
      }
    )
  })

  # Mostrar resultados de la simulación
  output$simulation_results <- renderDT({
    results <- simulation_results()
    req(!is.null(results), length(results) > 0)

    selected_sim <- get_selected_network_index(results)
    req(!is.null(selected_sim), selected_sim >= 1, selected_sim <= length(results))

    datatable(head(results[[selected_sim]], 1000),
      options = list(scrollX = TRUE, scrollY = "400px")
    )
  })

  # Función auxiliar para asegurar que el archivo tenga extensión .csv
  ensure_csv_extension <- function(filename) {
    if (tools::file_ext(filename) != "csv") {
      filename <- paste0(filename, ".csv")
    }
    return(filename)
  }

  # Función auxiliar para asegurar que el archivo tenga extensión .rds
  ensure_rds_extension <- function(filename) {
    if (tools::file_ext(filename) != "rds") {
      filename <- paste0(filename, ".rds")
    }
    return(filename)
  }

  # Función de verificación de NPEs
  verify_npes <- function(trials, NPEs) {
    missing_npes <- list()
    for (trial_type in names(trials)) {
      for (trial in trials[[trial_type]]) {
        configuration_step <- unlist(strsplit(trial, ","))
        npes_in_trial <- unique(configuration_step[seq(1, length(configuration_step) - 1, 2)])
        missing <- setdiff(npes_in_trial, NPEs$NPE)
        if (length(missing) > 0) {
          missing_npes[[trial_type]] <- c(missing_npes[[trial_type]], missing)
        }
      }
    }
    return(missing_npes)
  }

  # Guardar redes
  observeEvent(input$save_networks, {
    req(simulation_results())
    showModal(modalDialog(
      title = "Save networks",
      selectInput("save_option", "Save options:",
        choices = c("All networks in one file", "Individual networks")
      ),
      textInput("save_networks_filename", "File name (without extension):"),
      textInput("save_networks_path", "Save path (optional):", value = current_workdir()),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_save_networks", "Save")
      )
    ))
  })

  observeEvent(input$confirm_save_networks, {
    req(input$save_networks_filename)
    tryCatch(
      {
        results <- simulation_results()
        target_dir <- resolve_directory(input$save_networks_path)
        if (input$save_option == "All networks in one file") {
          combined_results <- do.call(rbind, lapply(seq_along(results), function(i) {
            cbind(Network = i, results[[i]])
          }))
          filename <- file.path(target_dir, paste0(input$save_networks_filename, ".csv"))
          write.csv(combined_results, filename, row.names = FALSE)
        } else {
          for (i in seq_along(results)) {
            filename <- file.path(target_dir, paste0(input$save_networks_filename, "_Network_", i, ".csv"))
            write.csv(results[[i]], filename, row.names = FALSE)
          }
        }
        removeModal()
        showNotification("Networks saved successfully", type = "message")
      },
      error = function(e) {
        showNotification(paste("Error saving the networks:", e$message), type = "error")
      }
    )
  })

  # Procesamiento de datos para gráficos
  processed_data <- reactiveVal(NULL)

  observeEvent(input$graficar, {
    if (!exists("datos_simulacion", envir = .GlobalEnv) || length(get("datos_simulacion", envir = .GlobalEnv)) == 0) {
      showNotification("There is no simulation data to plot. Please run a simulation first.", type = "error")
      return()
    }

    all_data <- get("datos_simulacion", envir = .GlobalEnv)
    selected_sim <- get_selected_network_index(all_data)
    if (is.null(selected_sim)) {
      showNotification("Please select a valid network first.", type = "error")
      return()
    }

    data <- all_data[[selected_sim]]
    if (is.null(data) || !is.data.frame(data) || nrow(data) == 0) {
      showNotification("Selected network has no valid data. Please run the simulation again.", type = "error")
      return()
    }

    # Asegurarse de que 'Phase' sea ordenada por contingencias/simulación
    data$Phase <- as.character(data$Phase)
    phase_order <- get_phase_order(data$Phase)
    data$Phase <- factor(data$Phase, levels = phase_order)
    data$Trial <- as.numeric(as.character(data$Trial))
    data$TimeStep <- as.numeric(as.character(data$TimeStep))

    # Identificar columnas que no son Phase, Trial, o TimeStep
    other_cols <- setdiff(names(data), c("Phase", "Trial", "TimeStep"))

    # Procesar datos para activaciones y pesos
    long_data <- data %>%
      pivot_longer(cols = all_of(other_cols), names_to = "Variable", values_to = "Value")

    # Separar activaciones y pesos
    activations <- long_data %>%
      filter(!grepl("-", Variable)) %>%
      rename(Unit = Variable, Activation = Value)

    weights <- long_data %>%
      filter(grepl("-", Variable)) %>%
      rename(Connection = Variable, Weight = Value)

    processed_data(list(activations = activations, weights = weights))

    # Obtener datos de NPEs
    npes_data <- npes()

    # Ordenar unidades
    ordered_units <- get_unit_order(unique(activations$Unit), npes_data)

    updateSelectizeInput(session, "activations_units",
      choices = ordered_units,
      selected = ordered_units[npes_data$Layer[match(ordered_units, npes_data$NPE)] == "PrimaryMotor"]
    )

    # Identificar unidades PrimaryMotor
    primary_motor_units <- ordered_units[npes_data$Layer[match(ordered_units, npes_data$NPE)] == "PrimaryMotor"]

    updateSelectizeInput(session, "aggregate_units",
      choices = ordered_units,
      selected = primary_motor_units
    )
    updateSelectInput(session, "aggregate_phase",
      choices = phase_order,
      selected = if (length(phase_order) > 0) phase_order[1] else character(0)
    )
    updateSelectInput(session, "selected_timestep",
      choices = unique(activations$TimeStep),
      selected = max(unique(activations$TimeStep)) - 1
    )
    updateSelectInput(session, "aggregate_timestep",
      choices = unique(activations$TimeStep),
      selected = max(unique(activations$TimeStep)) - 1
    )

    # Actualizar opciones de conexiones
    updateSelectizeInput(session, "weights_connections",
      choices = unique(weights$Connection),
      selected = unique(weights$Connection)[1:min(5, length(unique(weights$Connection)))]
    )
  })

  # Gráfico de Activaciones (Resultados individuales)
  output$activations_plot <- renderPlotly({
    req(processed_data(), input$activations_units, input$selected_timestep)
    data <- processed_data()$activations %>%
      filter(Unit %in% input$activations_units, TimeStep == input$selected_timestep)

    if (nrow(data) == 0) {
      return(plot_ly() %>% add_annotations(text = "No data for selected units/timestep.", showarrow = FALSE))
    }

    p <- ggplot(data, aes(x = Trial, y = Activation, color = Unit, group = Unit)) +
      geom_line(size = 0.7) +
      facet_wrap(~Phase, scales = "free_x") +
      theme_minimal() +
      labs(
        title = paste("Activations per trial (Timestep", input$selected_timestep, ")"),
        x = "Trial", y = "Activation"
      )

    ggplotly(p, tooltip = c("x", "y", "colour")) %>%
      layout(legend = list(orientation = "h", y = -0.2)) %>%
      config(displayModeBar = TRUE)
  })

  # Gráfico de Pesos de Conexiones (Resultados individuales)
  output$weights_plot <- renderPlotly({
    req(processed_data(), input$weights_connections, input$selected_timestep)
    data <- processed_data()$weights %>%
      filter(Connection %in% input$weights_connections, TimeStep == input$selected_timestep)

    if (nrow(data) == 0) {
      return(plot_ly() %>% add_annotations(text = "There is no data to display", showarrow = FALSE))
    }

    p <- ggplot(data, aes(x = Trial, y = Weight, color = Connection, group = Connection)) +
      geom_line(size = 0.8) +
      facet_wrap(~Phase, scales = "free_x") +
      theme_minimal() +
      labs(
        title = paste("Connection weights per trial (Timestep", input$selected_timestep, ")"),
        x = "Trial", y = "Weights"
      )

    ggplotly(p, tooltip = c("x", "y", "colour")) %>%
      layout(legend = list(orientation = "h", y = -0.2)) %>%
      config(displayModeBar = TRUE)
  })

  # Gráfico de Medidas Agregadas
  output$aggregate_plot <- renderPlotly({
    req(
      processed_data(), input$aggregate_phase, input$aggregate_units,
      input$aggregate_timestep, input$aggregate_measure, input$aggregate_error
    )

    data <- processed_data()$activations %>%
      filter(
        Phase == input$aggregate_phase,
        TimeStep == as.numeric(input$aggregate_timestep),
        Unit %in% input$aggregate_units
      )

    if (nrow(data) == 0) {
      return(plot_ly() %>% add_annotations(text = "No data to display", showarrow = FALSE))
    }

    aggregated_data <- data %>%
      group_by(Unit) %>%
      summarise(
        Mean = mean(Activation, na.rm = TRUE),
        Median = median(Activation, na.rm = TRUE),
        SE = sd(Activation, na.rm = TRUE) / sqrt(n()),
        SD = sd(Activation, na.rm = TRUE),
        .groups = "drop"
      )

    y_value <- ifelse(input$aggregate_measure == "mean", "Mean", "Median")
    error_value <- ifelse(input$aggregate_error == "se", "SE", "SD")

    p <- plot_ly() %>%
      add_trace(
        data = aggregated_data, x = ~Unit, y = as.formula(paste0("~", y_value)), type = "bar",
        marker = list(color = "rgba(158,202,225,0.6)", line = list(color = "rgb(8,48,107)", width = 1.5)),
        error_y = list(type = "data", array = aggregated_data[[error_value]], visible = TRUE),
        name = "Average"
      ) %>%
      layout(
        title = paste("Measure used:", input$aggregate_measure),
        xaxis = list(title = "Unit"),
        yaxis = list(title = "Average of activations", range = c(0, 1)),
        showlegend = TRUE
      )

    return(p)
  })

  # Procesamiento de datos para gráficos generales
  processed_data_general <- reactiveVal(NULL)

  observeEvent(input$graficar_general, {
    if (!exists("datos_simulacion", envir = .GlobalEnv) || length(get("datos_simulacion", envir = .GlobalEnv)) == 0) {
      showNotification("There is no simulation data to plot. Please run a simulation first.", type = "error")
      return()
    }

    all_data <- get("datos_simulacion", envir = .GlobalEnv)

    # Combinar todos los datos de todas las redes
    combined_data <- bind_rows(all_data, .id = "Network")

    # Asegurarse de que 'Phase' sea ordenada por contingencias/simulación
    combined_data$Phase <- as.character(combined_data$Phase)
    phase_order_general <- get_phase_order(combined_data$Phase)
    combined_data$Phase <- factor(combined_data$Phase, levels = phase_order_general)
    combined_data$Trial <- as.numeric(as.character(combined_data$Trial))
    combined_data$TimeStep <- as.numeric(as.character(combined_data$TimeStep))

    # Identificar columnas que no son Red, Phase, Trial, o TimeStep
    other_cols <- setdiff(names(combined_data), c("Network", "Phase", "Trial", "TimeStep"))

    # Procesar datos para activaciones y pesos
    long_data <- combined_data %>%
      pivot_longer(cols = all_of(other_cols), names_to = "Variable", values_to = "Value")

    # Separar activaciones y pesos
    activations <- long_data %>%
      filter(!grepl("-", Variable)) %>%
      rename(Unit = Variable, Activation = Value)

    weights <- long_data %>%
      filter(grepl("-", Variable)) %>%
      rename(Connection = Variable, Weight = Value)

    processed_data_general(list(activations = activations, weights = weights))

    # Obtener datos de NPEs
    npes_data <- npes()

    # Ordenar unidades
    ordered_units <- get_unit_order(unique(activations$Unit), npes_data)

    updateSelectizeInput(session, "activations_units_general",
      choices = ordered_units,
      selected = ordered_units[npes_data$Layer[match(ordered_units, npes_data$NPE)] == "PrimaryMotor"]
    )

    # Identificar unidades PrimaryMotor
    primary_motor_units <- ordered_units[npes_data$Layer[match(ordered_units, npes_data$NPE)] == "PrimaryMotor"]

    updateSelectizeInput(session, "aggregate_units_general",
      choices = ordered_units,
      selected = primary_motor_units
    )
    updateSelectInput(session, "aggregate_phase_general",
      choices = phase_order_general,
      selected = if (length(phase_order_general) > 0) phase_order_general[1] else character(0)
    )
    updateSelectInput(session, "selected_timestep_general",
      choices = unique(activations$TimeStep),
      selected = max(unique(activations$TimeStep)) - 1
    )
    updateSelectInput(session, "aggregate_timestep_general",
      choices = unique(activations$TimeStep),
      selected = max(unique(activations$TimeStep)) - 1
    )

    # Actualizar opciones de conexiones
    updateSelectizeInput(session, "weights_connections_general",
      choices = unique(weights$Connection),
      selected = unique(weights$Connection)[1:min(5, length(unique(weights$Connection)))]
    )
  })

  # Gráfico de Activaciones (Todas las Redes)
  output$activations_plot_general <- renderPlotly({
    req(processed_data_general(), input$activations_units_general, input$selected_timestep_general)
    data <- processed_data_general()$activations %>%
      filter(Unit %in% input$activations_units_general, TimeStep == input$selected_timestep_general)

    if (nrow(data) == 0) {
      return(plot_ly() %>% add_annotations(text = "No data for selected units/timestep.", showarrow = FALSE))
    }

    # Calcular el promedio de activación por ensayo y unidad
    avg_data <- data %>%
      group_by(Phase, Trial, Unit) %>%
      summarise(
        Avg_Activation = mean(Activation, na.rm = TRUE),
        SD_Activation = sd(Activation, na.rm = TRUE),
        .groups = "drop"
      )

    p <- ggplot(avg_data, aes(x = Trial, y = Avg_Activation, color = Unit, group = Unit)) +
      geom_line(size = 1.2) +
      facet_wrap(~Phase, scales = "free_x") +
      theme_minimal() +
      labs(
        title = paste("Average activations per trial (Timestep", input$selected_timestep_general, ")"),
        x = "Trial", y = "Average activation"
      )

    ggplotly(p, tooltip = c("x", "y", "colour")) %>%
      layout(legend = list(orientation = "h", y = -0.2)) %>%
      config(displayModeBar = TRUE)
  })

  # Gráfico de Pesos de Conexiones (Todas las Redes)
  output$weights_plot_general <- renderPlotly({
    req(processed_data_general(), input$weights_connections_general, input$selected_timestep_general)
    data <- processed_data_general()$weights %>%
      filter(Connection %in% input$weights_connections_general, TimeStep == input$selected_timestep_general)

    if (nrow(data) == 0) {
      return(plot_ly() %>% add_annotations(text = "No data to display", showarrow = FALSE))
    }

    # Calcular el promedio de peso por ensayo y conexión
    avg_data <- data %>%
      group_by(Phase, Trial, Connection) %>%
      summarise(
        Avg_Weight = mean(Weight, na.rm = TRUE),
        SD_Weight = sd(Weight, na.rm = TRUE),
        .groups = "drop"
      )

    p <- ggplot(avg_data, aes(x = Trial, y = Avg_Weight, color = Connection, group = Connection)) +
      geom_line(size = 1.2) +
      facet_wrap(~Phase, scales = "free_x") +
      theme_minimal() +
      labs(
        title = paste("Average connection weights per trial (Timestep", input$selected_timestep_general, ")"),
        x = "Trial", y = "Average weight"
      )

    ggplotly(p, tooltip = c("x", "y", "colour")) %>%
      layout(legend = list(orientation = "h", y = -0.2)) %>%
      config(displayModeBar = TRUE)
  })

  # Gráfico de Medidas Agregadas (Todas las Redes)
  output$aggregate_plot_general <- renderPlotly({
    req(
      processed_data_general(), input$aggregate_phase_general, input$aggregate_units_general,
      input$aggregate_timestep_general, input$aggregate_measure_general, input$aggregate_error_general
    )

    data <- processed_data_general()$activations %>%
      filter(
        Phase == input$aggregate_phase_general,
        TimeStep == as.numeric(input$aggregate_timestep_general),
        Unit %in% input$aggregate_units_general
      )

    if (nrow(data) == 0) {
      return(plot_ly() %>% add_annotations(text = "No data to display", showarrow = FALSE))
    }

    # Calcular estadísticas por unidad
    aggregated_data <- data %>%
      group_by(Unit) %>%
      summarise(
        Mean = mean(Activation, na.rm = TRUE),
        Median = median(Activation, na.rm = TRUE),
        SE = sd(Activation, na.rm = TRUE) / sqrt(n()),
        SD = sd(Activation, na.rm = TRUE),
        .groups = "drop"
      )

    # Calcular valores individuales por red
    individual_data <- data %>%
      group_by(Network, Unit) %>%
      summarise(
        Value = ifelse(input$aggregate_measure_general == "mean",
          mean(Activation, na.rm = TRUE),
          median(Activation, na.rm = TRUE)
        ),
        .groups = "drop"
      )

    y_value <- ifelse(input$aggregate_measure_general == "mean", "Mean", "Median")
    error_value <- ifelse(input$aggregate_error_general == "se", "SE", "SD")

    # Crear un vector de unidades únicas ordenadas
    unique_units <- unique(aggregated_data$Unit)

    p <- plot_ly() %>%
      add_trace(
        data = aggregated_data, x = ~ factor(Unit, levels = unique_units), y = as.formula(paste0("~", y_value)), type = "bar",
        marker = list(color = "rgba(158,202,225,0.6)", line = list(color = "rgb(8,48,107)", width = 1.5)),
        error_y = list(type = "data", array = aggregated_data[[error_value]], visible = TRUE),
        name = "Average"
      ) %>%
      add_trace(
        data = individual_data, x = ~ factor(Unit, levels = unique_units), y = ~Value, type = "scatter", mode = "markers",
        marker = list(color = "white", size = 8, line = list(color = "black", width = 1)),
        name = "Individual networks"
      ) %>%
      layout(
        title = paste("Measure used:", input$aggregate_measure_general),
        xaxis = list(
          title = "Type of NPE",
          type = "category",
          categoryorder = "array",
          categoryarray = unique_units
        ),
        yaxis = list(title = "Average of activations", range = c(0, 1)),
        showlegend = TRUE,
        legend = list(orientation = "h", y = -0.2),
        barmode = "overlay"
      )

    return(p)
  })

  # Guardar Ensayos
  observeEvent(input$save_trials, {
    req(input$trials_file_name != "")
    tryCatch(
      {
        filename <- ensure_rds_extension(input$trials_file_name)
        full_path <- file.path(resolve_directory(input$trials_file_path), filename)

        # Guardar la lista de ensayos directamente como un objeto RDS
        saveRDS(trials(), file = full_path)
        showNotification("Trials saved successfully", type = "message")
      },
      error = function(e) {
        showNotification(paste("Error saving the trials:", e$message), type = "error", duration = NULL)
      }
    )
  })

  # Importar Ensayos
  observeEvent(input$import_trials, {
    req(input$trials_file_name != "")
    tryCatch(
      {
        filename <- ensure_rds_extension(input$trials_file_name)
        full_path <- file.path(resolve_directory(input$trials_file_path), filename)
        if (!file.exists(full_path)) {
          stop(paste("File does not exist:", full_path))
        }

        # Leer el archivo RDS
        imported_trials <- readRDS(full_path)

        trials(imported_trials)
        assign("trials", imported_trials, envir = .GlobalEnv)
        updateSelectInput(session, "trial_types", choices = names(imported_trials))
        showNotification("Trials imported successfully", type = "message")
      },
      error = function(e) {
        showNotification(paste("Error importing trials:", e$message), type = "error", duration = NULL)
      }
    )
  })

  # Guardar Contingencias
  observeEvent(input$save_contingencies, {
    req(input$contingencies_file_name != "")
    tryCatch(
      {
        filename <- ensure_csv_extension(input$contingencies_file_name)
        full_path <- file.path(resolve_directory(input$contingencies_file_path), filename)

        # Convertir el vector de contingencias a un dataframe
        contingencies_df <- data.frame(Contingency = contingencies(), stringsAsFactors = FALSE)

        # Guardar como CSV
        write.csv(contingencies_df, file = full_path, row.names = FALSE)
        showNotification("Contingencies saved successfully", type = "message")
      },
      error = function(e) {
        showNotification(paste("Error saving the contingencies:", e$message), type = "error", duration = NULL)
      }
    )
  })

  # Importar Contingencias
  observeEvent(input$import_contingencies, {
    req(input$contingencies_file_name != "")
    tryCatch(
      {
        filename <- ensure_csv_extension(input$contingencies_file_name)
        full_path <- file.path(resolve_directory(input$contingencies_file_path), filename)
        if (!file.exists(full_path)) {
          stop(paste("File does not exist:", full_path))
        }

        # Leer el CSV
        imported_contingencies_df <- read.csv(full_path, stringsAsFactors = FALSE)

        # Convertir el dataframe de vuelta a un vector
        if ("Contingency" %in% names(imported_contingencies_df)) {
          imported_contingencies <- imported_contingencies_df$Contingency
        } else if ("Contingencia" %in% names(imported_contingencies_df)) {
          imported_contingencies <- imported_contingencies_df$Contingencia
        } else {
          stop("The contingencies file must include a 'Contingency' column.")
        }

        contingencies(imported_contingencies)
        assign("contingencies", imported_contingencies, envir = .GlobalEnv)
        showNotification("Contingencies imported successfully", type = "message")
      },
      error = function(e) {
        showNotification(paste("Error importing contingencies:", e$message), type = "error", duration = NULL)
      }
    )
  })

  # Botón de cierre
  observeEvent(input$close_app, {
    showModal(modalDialog(
      title = "Confirm closure",
      "Are you sure you want to close the program?",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_close", "Yes, close", class = "btn-danger")
      )
    ))
  })

  observeEvent(input$confirm_close, {
    stopApp()
  })

  # Help button observers
  observeEvent(input$help_npe_name, {
    showModal(modalDialog(
      title = "Help: NPU Name",
      "Name the network unit. Avoid apostrophes or quotes, and use a short ID. Example: S.1"
    ))
  })

  observeEvent(input$help_npe_type, {
    showModal(modalDialog(
      title = "Help: NPU Type",
      "Define whether your network unit is excitatory or inhibitory."
    ))
  })

  observeEvent(input$help_npe_layer, {
    showModal(modalDialog(
      title = "Help: Layer",
      "Type of unit you are going to create. The dopaminergic unit and US are already created by default."
    ))
  })

  observeEvent(input$reminder_us_d, {
    showTooltip("reminder_us_d", "Remember that the connection between US and D must be with a maximum weight of 1.", placement = "right", trigger = "hover")
  })

  observeEvent(input$help_npe_file_name, {
    showModal(modalDialog(
      title = "Help: NPUs File Name",
      "Register a name for your created units. This file will be saved with this name."
    ))
  })

  observeEvent(input$help_npe_file_path, {
    showModal(modalDialog(
      title = "Help: NPUs Directory Path",
      "Optional field. Leave it blank to use the current folder, or use the Home shortcuts to set all paths."
    ))
  })

  observeEvent(input$help_conn_pre, {
    showModal(modalDialog(
      title = "Help: Source NPU",
      "Select the source NPE to connect."
    ))
  })

  observeEvent(input$help_conn_post, {
    showModal(modalDialog(
      title = "Help: Target NPU",
      "Select the destination NPE to connect."
    ))
  })

  observeEvent(input$help_add_connection, {
    showModal(modalDialog(
      title = "Help: Add Connection",
      "Save the current connection definition to the connections table."
    ))
  })

  observeEvent(input$help_conn_file_name, {
    showModal(modalDialog(
      title = "Help: Connections File Name",
      "Register a name for your created connections. This file will be saved with this name."
    ))
  })

  observeEvent(input$help_conn_file_path, {
    showModal(modalDialog(
      title = "Help: Connections Directory Path",
      "Optional field. Leave it blank to use the current folder, or use the Home shortcuts to set all paths."
    ))
  })


  observeEvent(input$help_trial_name, {
    showModal(modalDialog(
      title = "Help: Trial Type Name",
      "Define the name of the trial type you want to create."
    ))
  })

  observeEvent(input$help_num_moments, {
    showModal(modalDialog(
      title = "Help: Number of Time Steps",
      "Shows how many time moments each of your S' (primary sensory) units will have."
    ))
  })

  observeEvent(input$help_phase_name, {
    showModal(modalDialog(
      title = "Help: Phase or Condition Name",
      "Record the name of your phase or condition, for example: Training."
    ))
  })

  observeEvent(input$help_presentation_mode, {
    showModal(modalDialog(
      title = "Help: Trial Presentation Mode",
      "Choose how trials are presented. Random is typically used for training, and In bulk for test phases."
    ))
  })

  observeEvent(input$help_trial_types, {
    showModal(modalDialog(
      title = "Help: Trial Types to Present",
      "Select which trial types will be used in that contingency. Example: for training choose X.1 and X.2; for testing choose XY."
    ))
  })

  observeEvent(input$help_iti_trial, {
    showModal(modalDialog(
      title = "Help: Add Existing ITI",
      "If you added an ITI you can select it here and decide how many time points it will have as a minimum and maximum."
    ))
  })

  observeEvent(input$help_num_simulations, {
    showModal(modalDialog(
      title = "Help: Number of Networks",
      "Choose the number of networks you want to simulate."
    ))
  })

  observeEvent(input$help_sim_file_name, {
    showModal(modalDialog(
      title = "Help: Simulation File Name",
      "Choose what name your file will have."
    ))
  })

  observeEvent(input$help_sim_file_path, {
    showModal(modalDialog(
      title = "Help: Simulation Directory Path",
      "Optional field. Leave it blank to use the current folder, or use the Home shortcuts to set all paths."
    ))
  })

  observeEvent(input$help_trials_file_name, {
    showModal(modalDialog(
      title = "Help: Trials File Name",
      "Set the file name used to save or load the trials object."
    ))
  })

  observeEvent(input$help_trials_file_path, {
    showModal(modalDialog(
      title = "Help: Trials Directory Path",
      "Optional field. Leave it blank to use the current folder, or use the Home shortcuts to set all paths."
    ))
  })

  observeEvent(input$help_contingencies_file_name, {
    showModal(modalDialog(
      title = "Help: Contingencies File Name",
      "Set the file name used to save or load the contingencies file."
    ))
  })

  observeEvent(input$help_contingencies_file_path, {
    showModal(modalDialog(
      title = "Help: Contingencies Directory Path",
      "Optional field. Leave it blank to use the current folder, or use the Home shortcuts to set all paths."
    ))
  })
}

# Ejecutar la aplicación Shiny
shinyApp(ui = ui, server = server)
