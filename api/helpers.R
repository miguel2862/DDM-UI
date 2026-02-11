# DDM Helper Functions

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
  for (i in 1:length(phases)) {
    current.phase <- trimws(unlist(strsplit(phases[i], ",")))
    trial.types <- trimws(unlist(strsplit(current.phase[3], "/")))
    if (any(!(trial.types %in% names(trials)))) stop("One of more trial names do not match trial names in phases")
  }
  return(list(valid = TRUE, message = "No errors found in the specification of the contingencies"))
}
