# DDM Phenomenon Templates
# Based on published experiments and verified simulation configurations

get_all_templates <- function() {
  list(
    extinction = list(
      id = "extinction",
      name = "Extinction",
      description = "A conditioned response decreases when the CS is repeatedly presented without the US. The model learns the CS-US association during training, then the response extinguishes when the US is omitted.",
      category = "Basic",
      available = TRUE,
      phases = "Training (100 trials) -> Extinction (100 trials)",
      npeCount = 7,
      connectionCount = 6
    ),
    alcala_2017 = list(
      id = "alcala_2017",
      name = "Autoshaped Impulsivity (Alcala, 2017)",
      description = "Simulation of autoshaped impulsive choice: Smaller-Sooner (SS) vs Larger-Later (LL) reinforcement. A context stimulus (S'2) with intermediate activation modulates responding to immediate (S'1) and delayed (S'3) stimuli. Based on Alcala (2017).",
      category = "Choice",
      available = TRUE,
      phases = "Training: SS/LL (100 each, Random, ITI) -> Test (25 trials)",
      npeCount = 13,
      connectionCount = 13
    ),
    burgos_donahoe_blocking = list(
      id = "burgos_donahoe_blocking",
      name = "Blocking (Burgos & Donahoe, 2016)",
      description = "Kamin blocking paradigm: Prior conditioning to stimulus A blocks conditioning to redundant stimulus X when presented in compound AX+. Demonstrates that prior predictive validity prevents learning about a new cue. Based on Burgos & Donahoe (2016).",
      category = "Compound",
      available = TRUE,
      phases = "A+ (100) -> AX+ (100) -> Test X (25)",
      npeCount = 13,
      connectionCount = 17
    ),
    burgos_donahoe_compound = list(
      id = "burgos_donahoe_compound",
      name = "Compound Training (Burgos & Donahoe, 2016)",
      description = "Compound stimulus training: AX+ compound is trained simultaneously. Both stimuli share the credit for predicting the US, allowing comparison with the blocking paradigm. Control condition for the blocking experiment. Based on Burgos & Donahoe (2016).",
      category = "Compound",
      available = TRUE,
      phases = "AX+ (100) -> Test X (25)",
      npeCount = 13,
      connectionCount = 17
    ),
    burgos_donahoe_successive = list(
      id = "burgos_donahoe_successive",
      name = "Successive Training (Burgos & Donahoe, 2016)",
      description = "Successive conditioning: Stimulus A is trained first, then stimulus X is trained separately. Both acquire associative strength independently, providing a baseline for comparison with blocking. Based on Burgos & Donahoe (2016).",
      category = "Compound",
      available = TRUE,
      phases = "A+ (100) -> X+ (100) -> Test A (25)",
      npeCount = 13,
      connectionCount = 17
    ),
    burgos_2000 = list(
      id = "burgos_2000",
      name = "Extinction & Reacquisition (Burgos, 2000)",
      description = "Acquisition, extinction, and reacquisition of a conditioned response. Demonstrates faster reacquisition after extinction, showing that extinction does not fully erase original learning. Uses a complex network with 3 sensory inputs and multiple motor outputs. Based on Burgos (2000).",
      category = "Basic",
      available = TRUE,
      phases = "Training I1 (300) -> Extinction I2 (300) -> Reacquisition I1 (300)",
      npeCount = 14,
      connectionCount = 32
    ),
    acquisition = list(
      id = "acquisition",
      name = "Acquisition",
      description = "The most fundamental learning phenomenon: a neutral CS paired with a biologically significant US gradually acquires the ability to elicit a conditioned response. The network learns to predict the US from the CS through increasing connection weights.",
      category = "Basic",
      available = TRUE,
      phases = "Training CS+US (100 trials)",
      npeCount = 7,
      connectionCount = 6
    ),
    latent_inhibition = list(
      id = "latent_inhibition",
      name = "Latent Inhibition",
      description = "Pre-exposure to a CS without consequence retards later conditioning when that CS is paired with a US. The prior non-reinforced experience reduces the associability of the CS, slowing subsequent learning.",
      category = "Basic",
      available = TRUE,
      phases = "Pre-exposure CS alone (100) -> Training CS+US (100)",
      npeCount = 7,
      connectionCount = 6
    )
  )
}


get_template_data <- function(template_id) {

  # -- Extinction ---------------------------------------------------------------
  if (template_id == "extinction") {
    return(list(
      id = "extinction",
      name = "Extinction",
      available = TRUE,
      npes = data.frame(
        NPE = c("US", "D", "S1", "S..1", "H1", "M..1", "M.1"),
        Type = rep("Excitatory", 7),
        Layer = c("US", "Dopaminergic", "PrimarySensory", "AssociativeSensory", "Hippocampal", "AssociativeMotor", "PrimaryMotor"),
        Activation = rep(0, 7),
        Temporal.Summation = rep(0.1, 7),
        Activation.Decay = rep(0.1, 7),
        mu = rep(0.2, 7),
        sigma = rep(0.15, 7),
        logisSigma = rep(0.1, 7),
        stringsAsFactors = FALSE
      ),
      connections = data.frame(
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
      ),
      hasITI = c(FALSE, FALSE)
    ))
  }

  # -- Alcala, Burgos & Aguayo-Mendoza (2017) - Autoshaped Impulsivity ---------
  if (template_id == "alcala_2017") {
    return(list(
      id = "alcala_2017",
      name = "Autoshaped Impulsivity (Alcala, 2017)",
      available = TRUE,
      npes = data.frame(
        NPE = c("US", "D", "S'1", "S'2", "S'3", "S''1", "S''2", "H1", "H2", "M''1", "M''2", "M'1", "M'2"),
        Type = rep("Excitatory", 13),
        Layer = c("US", "Dopaminergic",
                  "PrimarySensory", "PrimarySensory", "PrimarySensory",
                  "AssociativeSensory", "AssociativeSensory",
                  "Hippocampal", "Hippocampal",
                  "AssociativeMotor", "AssociativeMotor",
                  "PrimaryMotor", "PrimaryMotor"),
        Activation = rep(0, 13),
        Temporal.Summation = rep(0.1, 13),
        Activation.Decay = rep(0.1, 13),
        mu = rep(0.2, 13),
        sigma = rep(0.15, 13),
        logisSigma = rep(0.1, 13),
        stringsAsFactors = FALSE
      ),
      connections = data.frame(
        PreSinapticNPE = c("S'1", "S'2", "S'3", "S'2", "S''1", "S''1", "S''2", "S''2", "M''1", "M''1", "M''2", "M''2", "US"),
        PostSinapticNPE = c("S''1", "S''1", "S''2", "S''2", "H1", "M''1", "H2", "M''2", "D", "M'1", "D", "M'2", "D"),
        Weight = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1.0),
        alpha = rep(0.5, 13),
        beta = rep(0.1, 13),
        alpha_prime = rep(0.5, 13),
        beta_prime = rep(0.1, 13),
        stringsAsFactors = FALSE
      ),
      trials = list(
        SS = c(
          "US,0.00,S'1,1.00,S'2,0.65,S'3,0.00,True",
          "US,0.00,S'1,1.00,S'2,0.65,S'3,0.00,True",
          "US,0.70,S'1,1.00,S'2,0.65,S'3,0.00,True"
        ),
        LL = c(
          "US,0.00,S'1,0.00,S'2,0.65,S'3,1.00,True",
          "US,0.00,S'1,0.00,S'2,0.65,S'3,1.00,True",
          "US,0.00,S'1,0.00,S'2,0.65,S'3,1.00,True",
          "US,0.00,S'1,0.00,S'2,0.65,S'3,1.00,True",
          "US,0.00,S'1,0.00,S'2,0.65,S'3,1.00,True",
          "US,1.00,S'1,0.00,S'2,0.65,S'3,1.00,True"
        ),
        Prueba = c(
          "US,0.00,S'1,1.00,S'2,0.65,S'3,1.00,False",
          "US,0.00,S'1,1.00,S'2,0.65,S'3,1.00,False",
          "US,0.00,S'1,1.00,S'2,0.65,S'3,1.00,False",
          "US,0.00,S'1,1.00,S'2,0.65,S'3,1.00,False",
          "US,0.00,S'1,1.00,S'2,0.65,S'3,1.00,False"
        ),
        IEEn = c(
          "US,0.00,S'1,0.00,S'2,0.65,S'3,0.00,True"
        )
      ),
      contingencies = c(
        "Entrenamiento, Random, SS/LL, 100-100, True, 30, 30, IEEn",
        "Prueba, In bulk, Prueba, 25, False"
      ),
      hasITI = c(TRUE, FALSE)
    ))
  }

  # -- Burgos & Donahoe (2016) - Blocking --------------------------------------
  if (template_id == "burgos_donahoe_blocking") {
    npes_bur <- data.frame(
      NPE = c("US", "D", "A", "C", "X", "S..1", "S..2", "H1", "H2", "M..1", "M..2", "M.1", "M.2"),
      Type = rep("Excitatory", 13),
      Layer = c("US", "Dopaminergic",
                "PrimarySensory", "PrimarySensory", "PrimarySensory",
                "AssociativeSensory", "AssociativeSensory",
                "Hippocampal", "Hippocampal",
                "AssociativeMotor", "AssociativeMotor",
                "PrimaryMotor", "PrimaryMotor"),
      Activation = rep(0, 13),
      Temporal.Summation = rep(0.1, 13),
      Activation.Decay = rep(0.1, 13),
      mu = rep(0.2, 13),
      sigma = rep(0.15, 13),
      logisSigma = rep(0.1, 13),
      stringsAsFactors = FALSE
    )
    conns_bur <- data.frame(
      PreSinapticNPE = c("A", "C", "C", "X", "S..1", "S..1", "S..2", "S..2", "S..1", "S..2", "M..1", "M..2", "M..1", "M..2", "US", "US", "US"),
      PostSinapticNPE = c("S..1", "S..1", "S..2", "S..2", "H1", "M..1", "H2", "M..2", "M..2", "M..1", "D", "D", "M.1", "M.2", "D", "M.1", "M.2"),
      Weight = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1.0, 1.0, 1.0),
      alpha = rep(0.5, 17),
      beta = rep(0.12, 17),
      alpha_prime = rep(0.5, 17),
      beta_prime = rep(0.12, 17),
      stringsAsFactors = FALSE
    )
    return(list(
      id = "burgos_donahoe_blocking",
      name = "Blocking (Burgos & Donahoe, 2016)",
      available = TRUE,
      npes = npes_bur,
      connections = conns_bur,
      trials = list(
        "A+" = c(
          "US,0.00,A,1.00,C,0.90,X,0.00,True",
          "US,0.00,A,1.00,C,0.90,X,0.00,True",
          "US,0.00,A,1.00,C,0.90,X,0.00,True",
          "US,0.00,A,1.00,C,0.90,X,0.00,True",
          "US,1.00,A,1.00,C,0.90,X,0.00,True"
        ),
        "AX+" = c(
          "US,0.00,A,1.00,C,0.90,X,1.00,True",
          "US,0.00,A,1.00,C,0.90,X,1.00,True",
          "US,0.00,A,1.00,C,0.90,X,1.00,True",
          "US,0.00,A,1.00,C,0.90,X,1.00,True",
          "US,1.00,A,1.00,C,0.90,X,1.00,True"
        ),
        "X TST" = c(
          "US,0.00,A,0.00,C,0.90,X,1.00,False",
          "US,0.00,A,0.00,C,0.90,X,1.00,False",
          "US,0.00,A,0.00,C,0.90,X,1.00,False",
          "US,0.00,A,0.00,C,0.90,X,1.00,False",
          "US,0.00,A,0.00,C,0.90,X,1.00,False"
        )
      ),
      contingencies = c(
        "Entrenamiento, Random, A+, 100, False",
        "Bloqueo, Random, AX+, 100, False",
        "Test, In bulk, X TST, 25, False"
      ),
      hasITI = c(FALSE, FALSE, FALSE)
    ))
  }

  # -- Burgos & Donahoe (2016) - Compound Training -----------------------------
  if (template_id == "burgos_donahoe_compound") {
    npes_bur <- data.frame(
      NPE = c("US", "D", "A", "C", "X", "S..1", "S..2", "H1", "H2", "M..1", "M..2", "M.1", "M.2"),
      Type = rep("Excitatory", 13),
      Layer = c("US", "Dopaminergic",
                "PrimarySensory", "PrimarySensory", "PrimarySensory",
                "AssociativeSensory", "AssociativeSensory",
                "Hippocampal", "Hippocampal",
                "AssociativeMotor", "AssociativeMotor",
                "PrimaryMotor", "PrimaryMotor"),
      Activation = rep(0, 13),
      Temporal.Summation = rep(0.1, 13),
      Activation.Decay = rep(0.1, 13),
      mu = rep(0.2, 13),
      sigma = rep(0.15, 13),
      logisSigma = rep(0.1, 13),
      stringsAsFactors = FALSE
    )
    conns_bur <- data.frame(
      PreSinapticNPE = c("A", "C", "C", "X", "S..1", "S..1", "S..2", "S..2", "S..1", "S..2", "M..1", "M..2", "M..1", "M..2", "US", "US", "US"),
      PostSinapticNPE = c("S..1", "S..1", "S..2", "S..2", "H1", "M..1", "H2", "M..2", "M..2", "M..1", "D", "D", "M.1", "M.2", "D", "M.1", "M.2"),
      Weight = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 1.0, 1.0, 1.0),
      alpha = rep(0.5, 17),
      beta = rep(0.12, 17),
      alpha_prime = rep(0.5, 17),
      beta_prime = rep(0.12, 17),
      stringsAsFactors = FALSE
    )
    return(list(
      id = "burgos_donahoe_compound",
      name = "Compound Training (Burgos & Donahoe, 2016)",
      available = TRUE,
      npes = npes_bur,
      connections = conns_bur,
      trials = list(
        "AX+" = c(
          "US,0.00,A,1.00,C,0.90,X,1.00,True",
          "US,0.00,A,1.00,C,0.90,X,1.00,True",
          "US,0.00,A,1.00,C,0.90,X,1.00,True",
          "US,0.00,A,1.00,C,0.90,X,1.00,True",
          "US,1.00,A,1.00,C,0.90,X,1.00,True"
        ),
        "X TST" = c(
          "US,0.00,A,0.00,C,0.90,X,1.00,False",
          "US,0.00,A,0.00,C,0.90,X,1.00,False",
          "US,0.00,A,0.00,C,0.90,X,1.00,False",
          "US,0.00,A,0.00,C,0.90,X,1.00,False",
          "US,0.00,A,0.00,C,0.90,X,1.00,False"
        )
      ),
      contingencies = c(
        "EntrenamientoComp, Random, AX+, 100, False",
        "Prueba, In bulk, X TST, 25, False"
      ),
      hasITI = c(FALSE, FALSE)
    ))
  }

  # -- Burgos & Donahoe (2016) - Successive Training ---------------------------
  if (template_id == "burgos_donahoe_successive") {
    npes_bur <- data.frame(
      NPE = c("US", "D", "A", "C", "X", "S..1", "S..2", "H1", "H2", "M..1", "M..2", "M.1", "M.2"),
      Type = rep("Excitatory", 13),
      Layer = c("US", "Dopaminergic",
                "PrimarySensory", "PrimarySensory", "PrimarySensory",
                "AssociativeSensory", "AssociativeSensory",
                "Hippocampal", "Hippocampal",
                "AssociativeMotor", "AssociativeMotor",
                "PrimaryMotor", "PrimaryMotor"),
      Activation = rep(0, 13),
      Temporal.Summation = rep(0.1, 13),
      Activation.Decay = rep(0.1, 13),
      mu = rep(0.2, 13),
      sigma = rep(0.15, 13),
      logisSigma = rep(0.1, 13),
      stringsAsFactors = FALSE
    )
    conns_bur <- data.frame(
      PreSinapticNPE = c("A", "C", "C", "X", "S..1", "S..1", "S..2", "S..2", "S..1", "S..2", "M..1", "M..2", "M..1", "M..2", "US", "US", "US"),
      PostSinapticNPE = c("S..1", "S..1", "S..2", "S..2", "H1", "M..1", "H2", "M..2", "M..2", "M..1", "D", "D", "M.1", "M.2", "D", "M.1", "M.2"),
      Weight = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 1.0, 1.0, 1.0),
      alpha = rep(0.5, 17),
      beta = rep(0.12, 17),
      alpha_prime = rep(0.5, 17),
      beta_prime = rep(0.12, 17),
      stringsAsFactors = FALSE
    )
    return(list(
      id = "burgos_donahoe_successive",
      name = "Successive Training (Burgos & Donahoe, 2016)",
      available = TRUE,
      npes = npes_bur,
      connections = conns_bur,
      trials = list(
        "A+" = c(
          "US,0.00,A,1.00,C,0.90,X,0.00,True",
          "US,0.00,A,1.00,C,0.90,X,0.00,True",
          "US,0.00,A,1.00,C,0.90,X,0.00,True",
          "US,0.00,A,1.00,C,0.90,X,0.00,True",
          "US,1.00,A,1.00,C,0.90,X,0.00,True"
        ),
        "X+" = c(
          "US,0.00,A,0.00,C,0.90,X,1.00,True",
          "US,0.00,A,0.00,C,0.90,X,1.00,True",
          "US,0.00,A,0.00,C,0.90,X,1.00,True",
          "US,0.00,A,0.00,C,0.90,X,1.00,True",
          "US,1.00,A,0.00,C,0.90,X,1.00,True"
        ),
        "A TST" = c(
          "US,0.00,A,1.00,C,0.90,X,0.00,False",
          "US,0.00,A,1.00,C,0.90,X,0.00,False",
          "US,0.00,A,1.00,C,0.90,X,0.00,False",
          "US,0.00,A,1.00,C,0.90,X,0.00,False",
          "US,0.00,A,1.00,C,0.90,X,0.00,False"
        )
      ),
      contingencies = c(
        "Ent A+, Random, A+, 100, False",
        "Ent X+, Random, X+, 100, False",
        "A TST, In bulk, A TST, 25, False"
      ),
      hasITI = c(FALSE, FALSE, FALSE)
    ))
  }

  # -- Burgos (2000) - Extinction & Reacquisition ------------------------------
  if (template_id == "burgos_2000") {
    return(list(
      id = "burgos_2000",
      name = "Extinction & Reacquisition (Burgos, 2000)",
      available = TRUE,
      npes = data.frame(
        NPE = c("US", "D", "I1", "I2", "I3", "S..1", "S..2", "S..3", "H1", "M..1", "M..2", "M..3", "R1", "CR/UR"),
        Type = rep("Excitatory", 14),
        Layer = c("US", "Dopaminergic",
                  "PrimarySensory", "PrimarySensory", "PrimarySensory",
                  "AssociativeSensory", "AssociativeSensory", "AssociativeSensory",
                  "Hippocampal",
                  "AssociativeMotor", "AssociativeMotor", "AssociativeMotor",
                  "PrimaryMotor", "PrimaryMotor"),
        Activation = rep(0, 14),
        Temporal.Summation = rep(0.1, 14),
        Activation.Decay = rep(0.05, 14),
        mu = rep(0, 14),
        sigma = rep(1, 14),
        logisSigma = rep(0.1, 14),
        stringsAsFactors = FALSE
      ),
      connections = data.frame(
        PreSinapticNPE = c("I1", "I1", "I1", "I2", "I2", "I2",
                            "S..1", "S..1", "S..1", "S..1",
                            "S..2", "S..2", "S..2", "S..2",
                            "S..3", "S..3", "S..3", "S..3",
                            "M..1", "M..2", "M..3",
                            "M..1", "M..1", "M..2", "M..2", "M..3", "M..3",
                            "I3", "I3", "I3",
                            "US", "D"),
        PostSinapticNPE = c("S..1", "S..2", "S..3", "S..1", "S..2", "S..3",
                              "M..1", "M..2", "M..3", "H1",
                              "H1", "M..1", "M..2", "M..3",
                              "H1", "M..1", "M..2", "M..3",
                              "D", "D", "D",
                              "R1", "CR/UR", "R1", "CR/UR", "R1", "CR/UR",
                              "S..1", "S..2", "S..3",
                              "D", "CR/UR"),
        Weight = c(rep(0.01, 30), 1.0, 1.0),
        alpha = rep(0.5, 32),
        beta = rep(0.035, 32),
        alpha_prime = rep(0.035, 32),
        beta_prime = rep(0.1, 32),
        stringsAsFactors = FALSE
      ),
      trials = list(
        I1 = c(
          "US,0.00,I1,1.00,I2,0.00,I3,0.00,True",
          "US,0.00,I1,1.00,I2,0.00,I3,0.00,True",
          "US,0.00,I1,1.00,I2,0.00,I3,0.00,True",
          "US,0.00,I1,1.00,I2,0.00,I3,0.00,True",
          "US,0.00,I1,1.00,I2,0.00,I3,0.00,True",
          "US,1.00,I1,1.00,I2,0.00,I3,0.00,True"
        ),
        I2 = c(
          "US,0.00,I1,0.00,I2,1.00,I3,0.00,True",
          "US,0.00,I1,0.00,I2,1.00,I3,0.00,True",
          "US,0.00,I1,0.00,I2,1.00,I3,0.00,True",
          "US,0.00,I1,0.00,I2,1.00,I3,0.00,True",
          "US,0.00,I1,0.00,I2,1.00,I3,0.00,True",
          "US,0.00,I1,0.00,I2,1.00,I3,0.00,True"
        )
      ),
      contingencies = c(
        "Entrenamiento, Random, I1, 300, False",
        "Extincion, Random, I2, 300, False",
        "Readquisicion, Random, I1, 300, False"
      ),
      hasITI = c(FALSE, FALSE, FALSE)
    ))
  }

  # -- Acquisition (basic CS-US pairing) -----------------------------------------
  if (template_id == "acquisition") {
    return(list(
      id = "acquisition",
      name = "Acquisition",
      available = TRUE,
      npes = data.frame(
        NPE = c("US", "D", "S1", "S..1", "H1", "M..1", "M.1"),
        Type = rep("Excitatory", 7),
        Layer = c("US", "Dopaminergic", "PrimarySensory", "AssociativeSensory", "Hippocampal", "AssociativeMotor", "PrimaryMotor"),
        Activation = rep(0, 7),
        Temporal.Summation = rep(0.1, 7),
        Activation.Decay = rep(0.1, 7),
        mu = rep(0.2, 7),
        sigma = rep(0.15, 7),
        logisSigma = rep(0.1, 7),
        stringsAsFactors = FALSE
      ),
      connections = data.frame(
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
        )
      ),
      contingencies = c(
        "training, Random, Training, 100, False"
      ),
      hasITI = c(FALSE)
    ))
  }


  # -- Latent Inhibition (pre-exposure CS alone -> acquisition CS+US) ----------
  if (template_id == "latent_inhibition") {
    return(list(
      id = "latent_inhibition",
      name = "Latent Inhibition",
      available = TRUE,
      npes = data.frame(
        NPE = c("US", "D", "S1", "S..1", "H1", "M..1", "M.1"),
        Type = rep("Excitatory", 7),
        Layer = c("US", "Dopaminergic", "PrimarySensory", "AssociativeSensory", "Hippocampal", "AssociativeMotor", "PrimaryMotor"),
        Activation = rep(0, 7),
        Temporal.Summation = rep(0.1, 7),
        Activation.Decay = rep(0.1, 7),
        mu = rep(0.2, 7),
        sigma = rep(0.15, 7),
        logisSigma = rep(0.1, 7),
        stringsAsFactors = FALSE
      ),
      connections = data.frame(
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
        PreExposure = c(
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True"
        ),
        Training = c(
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True",
          "US,0.00,S1,1.00,True",
          "US,1.00,S1,1.00,True"
        )
      ),
      contingencies = c(
        "pre-exposure, Random, PreExposure, 100, False",
        "training, Random, Training, 100, False"
      ),
      hasITI = c(FALSE, FALSE)
    ))
  }

  return(list(error = "Template not found or not yet available"))
}
