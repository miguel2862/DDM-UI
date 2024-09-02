#################################################################################################
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
# Cargar las funciones del modelo DDM
Simulate.DBP <- function(NPEs, Connections, TimeSteps, HasITI){
  
  tryCatch({
    
    setClass("NPE", slots = c(Activation = "numeric",
                              PreviousActivation = "numeric",
                              ActivationDecay = "numeric",
                              ExcitatoryInput = "numeric",
                              InhibitoryInput="numeric",
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
                              r = "numeric"))
    
    setClass("Connection", slots = c(weight = "numeric",
                                     Name = "character",
                                     alpha = "numeric",
                                     beta = "numeric",
                                     alpha_prime = "numeric",
                                     beta_prime = "numeric",
                                     preSinapticNPE = "character",
                                     p = "numeric"))
    
    
    L <- function(x, sigma){
      return(1 / (1 + exp((-x + 0.5) / sigma)))
    }
    
    dVTA <- function(){
      d = 0
      n = 0
      for(i in 1:length(network)){
        if(network[[i]]@Layer == "Dopaminergic"){
          d = d + (network[[i]]@Activation - network[[i]]@PreviousActivation)
          n = n + 1
        }
      }
      return(ifelse(d>0,d/n,0))
    }
    
    dCA1 <- function(dVTA){
      d = 0
      n = 0
      for(i in 1:length(network)){
        if(network[[i]]@Layer == "Hippocampal"){
          d = d + abs(network[[i]]@Activation-network[[i]]@PreviousActivation)
          #d = d + (network[[i]]@Activation-network[[i]]@PreviousActivation)
          n = n + 1
        }
      }
      return(d/n+dVTA*(1-d/n))
    }
    
    Compute.r <- function(npe){
      sum.weights.exc <- 0
      sum.weights.inh <- 0
      
      if(length(network[[npe]]@InputConnections)>0){
        for(i in 1:length(network[[npe]]@InputConnections)){
          pre <- network[[npe]]@InputConnections[[i]]@preSinapticNPE
          
          if (network[[pre]]@Type == "Excitatory"){
            sum.weights.exc <- sum.weights.exc + ifelse(network[[pre]]@Layer == "US",0,network[[npe]]@InputConnections[[i]]@weight)
          }else{
            sum.weights.inh <- sum.weights.inh + network[[npe]]@InputConnections[[i]]@weight
          }
        }
      }
      
      return(c(1-sum.weights.exc, 1-sum.weights.inh))
    }
    
    data.sim <- as.data.frame(matrix(nrow = nrow(TimeSteps), ncol = nrow(NPEs) + nrow(Connections) + 3))
    
    colnames(data.sim) <- c("Phase","Trial","TimeStep",NPEs[,1],paste(Connections[,1],Connections[,2], sep = "-"))
    
    ### Initialize network with NPEs and connections
    network <- list()
    
    for(i in 1:nrow(NPEs)){
      network[[NPEs[i,1]]] <- new("NPE", Name=NPEs[i,1],
                                  Type = NPEs[i,2],
                                  Layer = NPEs[i,3],
                                  Activation = NPEs[i,4],
                                  TemporalSummation = NPEs[i,5],
                                  ActivationDecay = NPEs[i,6],
                                  mu = NPEs[i,7],
                                  sigma = NPEs[i,8],
                                  logisSigma = NPEs[i,9],
                                  PreviousExcitatoryInput = 0,
                                  ExcitatoryInput = 0,
                                  InhibitoryInput = 0,
                                  PreviousActivation = 0)
    }
    
    for(i in 1:nrow(Connections)){
      connection.name <- paste(Connections[i,1],Connections[i,2],sep = "-")
      currentPostSinapticNPE <- Connections[i,2]
      network[[currentPostSinapticNPE]]@InputConnections[[connection.name]] <- new("Connection", Name=connection.name,
                                                                                   weight = Connections[i,3],
                                                                                   alpha = Connections[i,4],
                                                                                   beta = 0.12,  # Valor por defecto de beta
                                                                                   alpha_prime = Connections[i,6],
                                                                                   beta_prime = 0.12,  # Valor por defecto de beta prime
                                                                                   preSinapticNPE = Connections[i,1])
    }
    
    LearningRuleIsActive = F
    PreviousdCA1 = 0
    
    pb <- txtProgressBar(max = nrow(TimeSteps), style = 3)
    
    current.phase <- 1
    
    for(ts in 1:nrow(TimeSteps)){
      
      ### save phase, trial and timestep
      data.sim[ts,1] <- TimeSteps[ts,1]
      
      if(ts>1 && TimeSteps[ts-1,1]!=TimeSteps[ts,1]){
        current.phase = current.phase + 1
      }
      
      data.sim[ts,2] <- TimeSteps[ts,2]
      data.sim[ts,3] <- TimeSteps[ts,3]
      
      # reset activations if HasITI[current.phase] is false and trial onset
      if(TimeSteps[ts,3]==1 & !HasITI[current.phase]){
        for(i in 1:length(network)){
          network[[i]]@Activation <- L(0, network[[i]]@logisSigma)
          network[[i]]@ExcitatoryInput <- L(0, network[[i]]@logisSigma)
        }
      }
      
      ### reset input activations to zero
      for (npu in 1:length(network)){
        if(network[[npu]]@Layer == 'US' | network[[npu]]@Layer == 'PrimarySensory') {
          network[[npu]]@Activation = 0
        }
      } 
      
      # set inputs
      LearningRuleIsActive <- as.logical(TimeSteps[ts,ncol(TimeSteps)])
      
      for(unit in seq(4,ncol(TimeSteps)-1,2)){
        network[[TimeSteps[ts,unit]]]@Activation <- TimeSteps[ts,unit+1]
      }
      
      # update activations to all members of the network list (scrambled) # asynchronous random
      scrambledNPEs <- sample(1:length(network),length(network),replace = F)
      
      for(i in scrambledNPEs){
        currentPostSinapticNPE <- network[[i]]@Name
        
        network[[i]]@PreviousActivation <- network[[i]]@Activation
        network[[i]]@PreviousExcitatoryInput = network[[i]]@ExcitatoryInput
        
        if(network[[i]]@Layer=="US" | network[[i]]@Layer=="PrimarySensory") {
          data.sim[ts,network[[i]]@Name] <- network[[i]]@Activation
          next
        }
        
        #checks whether current NPE is connected to a US unit
        npe.is.connected.to.an.active.US <- F
        
        if(length(network[[i]]@InputConnections)>0){
          for(j in 1:length(network[[i]]@InputConnections)){
            pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
            if(network[[pre]]@Layer == "US" && network[[pre]]@Activation >0){
              npe.is.connected.to.an.active.US <- T
              break
            }
          }
        }
        
        if(npe.is.connected.to.an.active.US){
          # unconditional activation
          network[[i]]@Activation = network[[pre]]@Activation
        }
        else{
          # conditional activation
          network[[i]]@Threshold <- rnorm(1,network[[i]]@mu,network[[i]]@sigma)
          
          # compute excitatory and inhibitory inputs
          if(length(network[[i]]@InputConnections)>0){
            value <- c(0,0)
            
            for(j in 1:length(network[[i]]@InputConnections)){
              pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
              if(network[[pre]]@Type == "Excitatory"){
                value[1] = value[1] + network[[pre]]@Activation * network[[i]]@InputConnections[[j]]@weight
              }else{
                value[2] = value[2] + network[[pre]]@Activation * network[[i]]@InputConnections[[j]]@weight
              }
            }
            network[[i]]@ExcitatoryInput = value[1]
            network[[i]]@InhibitoryInput = value[2]
          }else{
            network[[i]]@ExcitatoryInput = 0
            network[[i]]@InhibitoryInput = 0
          }
          
          p_epsp = L(network[[i]]@ExcitatoryInput, network[[i]]@logisSigma)
          
          if(is.nan(p_epsp)) browser()
          
          p_ipsp = ifelse(network[[i]]@InhibitoryInput == 0, 0, L(network[[i]]@InhibitoryInput, network[[i]]@logisSigma))
          
          if(network[[i]]@Layer!= "PrimarySensory" & network[[i]]@Layer!= "US"){
            if(p_epsp > p_ipsp){
              if(p_epsp >= network[[i]]@Threshold){
                #Reactivation
                network[[i]]@Activation = p_epsp + network[[i]]@TemporalSummation * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) * (1 - p_epsp) - p_ipsp
              } else{
                # decay
                network[[i]]@Activation = L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma) - network[[i]]@ActivationDecay * L(network[[i]]@PreviousExcitatoryInput, network[[i]]@logisSigma)
              }
            }else{
              #inhibition
              network[[i]]@Activation = 0
            }
          }else{
            network[[i]]@Activation = network[[i]]@Activation
          }
        }
        
        data.sim[ts,network[[i]]@Name] <- network[[i]]@Activation
      }
      
      # update weights to all members of the inputConnections list of all the NPEs
      if(LearningRuleIsActive){
        dD = dVTA()
        dH = dCA1(dD)
        PreviousdCA1 = dH
        
        scrambledNPEs <- sample(1:length(network),length(network),replace = F)
        
        for(i in scrambledNPEs){
          currentPostSinapticNPE <- network[[i]]@Name
          
          if(network[[i]]@Layer=="US" | network[[i]]@Layer=="PrimarySensory") {next}
          
          network[[i]]@r <- Compute.r(i)
          
          scrambleConnections <- sample(1:length(network[[i]]@InputConnections),length(network[[i]]@InputConnections),replace = F)
          
          if(network[[i]]@Layer %in% c("AssociativeSensory","Hippocampal")){
            d = dH
          }else{
            d = dD
          }
          
          for (j in scrambleConnections){
            pre <- network[[i]]@InputConnections[[j]]@preSinapticNPE
            
            if(network[[pre]]@Layer != "US"){
              if(d>0.001){
                if(network[[pre]]@Type == 'Excitatory'){
                  network[[i]]@InputConnections[[j]]@p <- ifelse(network[[i]]@ExcitatoryInput == 0, 0, network[[pre]]@Activation * network[[i]]@InputConnections[[j]]@weight / network[[i]]@ExcitatoryInput)
                  
                  network[[i]]@InputConnections[[j]]@weight = network[[i]]@InputConnections[[j]]@weight + 
                    network[[i]]@InputConnections[[j]]@alpha*network[[i]]@r[1]*
                    network[[i]]@Activation*d*network[[i]]@InputConnections[[j]]@p
                }else{
                  network[[i]]@InputConnections[[j]]@p <- ifelse(network[[i]]@InhibitoryInput == 0, 0, network[[pre]]@Activation * network[[i]]@InputConnections[[j]]@weight / network[[i]]@InhibitoryInput)
                  
                  network[[i]]@InputConnections[[j]]@weight = network[[i]]@InputConnections[[j]]@weight + 
                    network[[i]]@InputConnections[[j]]@alpha_prime*network[[i]]@r[2]*
                    network[[i]]@Activation*d*network[[i]]@InputConnections[[j]]@p
                }
              }
              else{
                network[[i]]@InputConnections[[j]]@weight = network[[i]]@InputConnections[[j]]@weight - 
                  ifelse(network[[pre]]@Type == "Excitatory",
                         network[[i]]@InputConnections[[j]]@beta, 
                         network[[i]]@InputConnections[[j]]@beta_prime) *
                  network[[i]]@InputConnections[[j]]@weight * network[[pre]]@Activation * 
                  network[[i]]@Activation
              }
            }
            
            data.sim[ts,network[[i]]@InputConnections[[j]]@Name] <- network[[i]]@InputConnections[[j]]@weight
          }
        }
      }else{
        for(i in scrambledNPEs){
          scrambleConnections <- sample(1:length(network[[i]]@InputConnections),length(network[[i]]@InputConnections),replace = F)
          
          for (j in scrambleConnections){
            data.sim[ts,network[[i]]@InputConnections[[j]]@Name] <- network[[i]]@InputConnections[[j]]@weight
          }
        }
      }
      
      setTxtProgressBar(pb, ts) 
    }
    return(data.sim)
  }, 
  error = function(e) e)
}
##############################################################################################
Create.Phases <- function(phases, trials){
  
  # phases is a character vector with comma delimited characters. For example,
  # "Training 1, Random,A+/X-,100-100,False,30,30,ITI entrenamiento"
  # "Training 2, In bulk,AX+,100,False,30,30,ITI entrenamiento"
  # "Test,       In bulk,Prueba X,25,False,30,30,ITI entrenamiento"
  
  #trials is a list with named elements
  # each element is a named vector of comma delimited characters
  # for example, trials[["A+"]] might be:
  
  #   [1] "S.1,1,S.2,0.65,S.3,0,US,0,True"
  #   [2] "S.1,1,S.2,0.65,S.3,0,US,0,True"
  #   [3] "S.1,1,S.2,0.65,S.3,0,US,0,True"
  #   [4] "S.1,1,S.2,0.65,S.3,0,US,0,True"
  #   [5] "S.1,1,S.2,0.65,S.3,0,US,1,True"
  
  timesteps <- as.data.frame(matrix(nrow=0,ncol = 3 + length(unlist(strsplit(trials[[1]][1],",")))))
  
  for(i in 1:length(phases)){
    
    current.ts <- 1
    current.trial <- 1
    
    current.phase <- trimws(unlist(strsplit(phases[i],",")))
    
    phase.name <- current.phase[1]
    trial.order <- tolower(current.phase[2])
    trial.types <- trimws(unlist(strsplit(current.phase[3],"/")))
    
    if(any(!(trial.types %in% names(trials)))) stop("One of more trial names do not match trial names in phases")
    
    trial.numbers <- as.integer(trimws(unlist(strsplit(current.phase[4],"-"))))
    has.iti <- as.logical(trimws(current.phase[5]))
    
    if(has.iti){
      min.ITI <- as.integer(trimws(current.phase[6]))
      max.ITI <- as.integer(trimws(current.phase[7]))
      ITITimestep <- trimws(current.phase[8])
      
    }
    
    if(trial.order == "in bulk"){
      
      for(j in 1:length(trial.types)){
        
        for(k in 1:trial.numbers[j]){
          
          if(has.iti){
            
            if(!(ITITimestep %in% names(trials))){ stop("ITI time steps not in trials")}
            
            current.ITI <- ifelse(min.ITI == max.ITI, min.ITI, sample(min.ITI:max.ITI, 1))
            
            for (ts in 1:current.ITI){
              
              timesteps<- rbind(timesteps,c(phase.name,current.trial, current.ts,trimws(unlist(strsplit(trials[[ITITimestep]],",")))))
              current.ts=current.ts+1
            }
            
          }else{current.ts <- 1}
          
          for(ts in 1:length(trials[[trial.types[j]]])){
            timesteps<- rbind(timesteps,c(phase.name,current.trial, current.ts, trimws(unlist(strsplit(trials[[trial.types[j]]][ts],",")))))
            current.ts=current.ts+1
          }
          current.trial<- current.trial+1
          current.ts <- 1
        }
        
      }
      
      
    }else if(trial.order == "alternated"){
      
      for(j in 1:trial.numbers[1]){
        
        for(k in 1:length(trial.types)){
          
          if(has.iti){
            
            if(!(ITITimestep %in% names(trials))){ stop("ITI time steps not in trials")}
            
            current.ITI <- ifelse(min.ITI == max.ITI, min.ITI, sample(min.ITI:max.ITI, 1))
            
            for (ts in 1:current.ITI){
              
              timesteps<- rbind(timesteps,c(phase.name,current.trial, current.ts,trimws(unlist(strsplit(trials[[ITITimestep]],",")))))
              current.ts=current.ts+1
            }
            
          }else{current.ts <- 1}
          
          for(ts in 1:length(trials[[trial.types[k]]])){
            timesteps<- rbind(timesteps,c(phase.name,current.trial, current.ts, trimws(unlist(strsplit(trials[[trial.types[k]]][ts],",")))))
            current.ts=current.ts+1
          }
          current.trial<- current.trial+1
          current.ts <- 1
        }
        
      }
      
    }else if(trial.order == "random"){
      
      trial.sequence <- vector(mode = "integer")
      
      for(n in 1:length(trial.numbers)){
        trial.sequence <- c(trial.sequence,rep(n,trial.numbers[n]))
      }
      
      trial.sequence <- sample(trial.sequence,length(trial.sequence),replace = F)
      
      for(j in trial.sequence){
        
        if(has.iti){
          
          if(!(ITITimestep %in% names(trials))){ stop("ITI time steps not in trials")}
          
          current.ITI <- ifelse(min.ITI == max.ITI, min.ITI, sample(min.ITI:max.ITI, 1))
          
          for (ts in 1:current.ITI){
            
            timesteps<- rbind(timesteps,c(phase.name,current.trial, current.ts,trimws(unlist(strsplit(trials[[ITITimestep]],",")))))
            current.ts=current.ts+1
          }
          
        }else{current.ts <- 1}
        
        for(ts in 1:length(trials[[trial.types[j]]])){
          timesteps<- rbind(timesteps,c(phase.name,current.trial,current.ts, trimws(unlist(strsplit(trials[[trial.types[j]]][ts],",")))))
          current.ts=current.ts+1
        }
        current.trial<- current.trial+1
        current.ts <- 1
        
      }
      
    }
    
    
  }
  
  columns <- suppressWarnings(which(!is.na(sapply(timesteps[1,], as.numeric))==T))
  
  timesteps[,columns]<-as.data.frame(apply(timesteps[,columns],2,as.numeric))
  
  print("TimeSteps created")
  
  return(timesteps)
  
  
  
}

create.NPEs <- function(npe){
  
  NPEs <- as.data.frame(matrix(nrow = length(npe), ncol = 9))
  
  for(n in 1:length(npe)){
    NPEs[n,] <- trimws(unlist(strsplit(npe[n],",")))
  }
  columns <- suppressWarnings(which(!is.na(sapply(NPEs[1,], as.numeric))==T))
  
  NPEs[,columns]<-as.data.frame(apply(NPEs[,columns],2,as.numeric))
  colnames(NPEs) <- c("NPE", "Type", "Layer","Activation", "Temporal.Summation","Activation.Decay", "mu", "sigma", "logisSigma" )
  return(NPEs)
}

create.Connections <- function(conn, NPEs){
  
  Connections <- as.data.frame(matrix(nrow = length(conn), ncol = 7))
  
  for(n in 1:length(conn)){
    Connections[n,] <- trimws(unlist(strsplit(conn[n],",")))
    
    if(any(!Connections[n,1:2] %in% NPEs)) stop("Either the preSinaptic or the postSinaptic NPE is not a member of the NPEs")
    
  }
  columns <- suppressWarnings(which(!is.na(sapply(Connections[1,], as.numeric))==T))
  
  Connections[,columns]<-as.data.frame(apply(Connections[,columns],2,as.numeric))
  colnames(Connections) <- c("PreSinapticNPE", "PostSinapticNPE", "Weight", "alpha", "beta", "alpha_prime", "beta_prime")
  return(Connections)
}

#################################################################################################

ui <- function(request) {
  tagList(
    useShinyjs(),
    tags$script('$(function () { $("[data-toggle=\'tooltip\']").tooltip(); })'),
    tags$head(
      tags$style(HTML("
    /* Styles for the splash screen */
    #splash-screen {
        position: fixed;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background: linear-gradient(135deg, #1a237e 0%, #4a148c 100%);
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
        background-color: #f4f6f9;
    }
    .box {
        border-top: 3px solid #3c8dbc;
    }
    .box-title {
        font-weight: bold;
    }
    /* Style for the 'Create architecture' button */
    #create_architecture {
        background-color: #f39c12;
        color: white;
        font-weight: bold;
        font-size: 16px;
        padding: 10px 20px;
        margin-top: 20px;
        display: block;
        width: 100%;
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
    div(id = "splash-screen",
        div(id = "splash-content",
            tags$svg(id = "splash-logo", viewBox = "0 0 100 100",
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
            p(class = "splash-info", "Last modified: August 2024"),
            p(class = "splash-info", "University of Guadalajara")
        )
    ),
    # Main content of the application
    div(id = "main-content",
        dashboardPage(
          dashboardHeader(title = "DDM Simulator"),
          dashboardSidebar(
            sidebarMenu(
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
              style = "position: absolute; bottom: 0; left: 0; right: 0.4; padding: 10px;",
              actionButton("close_app", "Close Program", 
                           icon = icon("power-off"), 
                           style = "width: 100%; color: #fff; background-color: #d9534f; border-color: #d43f3a;")
            )
          ),
          dashboardBody(
            tabItems(
              # Home tab
              tabItem(tabName = "home",
                      fluidRow(
                        box(
                          title = "Welcome to the Diffuse Discrepancy Model (DDM) Simulator",
                          width = 12,
                          p("The Diffuse Discrepancy Model (DDM) is a powerful tool for simulating learning and conditioning phenomena in behavioral sciences."),
                          p("Originally developed by Donahoe, Burgos and Palmer (1993), the DDM offers a connectionist interpretation of the unified principle of reinforcement for operant and Pavlovian conditioning."),
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
                          p("For more information on how to use each component of the simulator, please refer to the 'Help' section.")
                        )
                      )
              ),
              # Network architecture tab
              tabItem(tabName = "network",
                      fluidRow(
                        box(
                          title = "Units (NPUs)",
                          width = 6,
                          div(style = "display: flex; align-items: center;",
                              textInput("npe_name", "NPU Name"),
                              icon("question-circle", id = "help_npe_name", class = "help-icon")
                          ),
                          fluidRow(
                            column(6, 
                                   div(style = "display: flex; align-items: center;",
                                       selectInput("npe_type", "NPU Type", 
                                                   choices = c("Excitatory", "Inhibitory"),
                                                   width = "100%"),
                                       icon("question-circle", id = "help_npe_type", class = "help-icon")
                                   )
                            ),
                            column(6,
                                   div(style = "display: flex; align-items: center;",
                                       selectInput("npe_layer", "Layer", 
                                                   choices = c("PrimarySensory", "AssociativeSensory", "Hippocampal", "AssociativeMotor", "PrimaryMotor", "Dopaminergic"),
                                                   width = "100%"),
                                       icon("question-circle", id = "help_npe_layer", class = "help-icon")
                                   )
                            )
                          ),
                          fluidRow(
                            column(4, numericInput("npe_activation", "Initial Activation", value = 0, min = 0, max = 1, step = 0.1)),
                            column(4, numericInput("npe_temporal_summation", "Temporal Summation (τ)", value = 0.1, min = 0, max = 1, step = 0.1)),
                            column(4, numericInput("npe_decay", "Decay Rate (κ)", value = 0.1, min = 0, max = 1, step = 0.1))
                          ),
                          fluidRow(
                            column(4, numericInput("npe_mu", "Threshold Mean (μ)", value = 0.2, min = 0, max = 1, step = 0.1)),
                            column(4, numericInput("npe_sigma", "Threshold Deviation (σ)", value = 0.15, min = 0, max = 1, step = 0.1)),
                            column(4, numericInput("npe_logistic_slope", "Logistic Slope", value = 0.1, min = 0, max = 1, step = 0.1))
                          ),
                          actionButton("add_npe", "Add NPU"),
                          hr(),
                          div(style = "display: flex; align-items: center;",
                              textInput("npe_file_name", "NPUs File Name"),
                              icon("question-circle", id = "help_npe_file_name", class = "help-icon")
                          ),
                          div(style = "display: flex; align-items: center;",
                              textInput("npe_file_path", "NPUs Directory Path"),
                              icon("question-circle", id = "help_npe_file_path", class = "help-icon")
                          ),
                          actionButton("save_npes", "Save NPUs"),
                          actionButton("import_npes", "Import NPUs")
                        ),
                        
                        box(
                          title = "Connections",
                          width = 6,
                          div(style = "display: flex; align-items: center;",
                              selectInput("conn_pre", "Source NPU", choices = NULL),
                              icon("question-circle", id = "help_conn_pre", class = "help-icon")
                          ),
                          div(style = "display: flex; align-items: center;",
                              selectInput("conn_post", "Target NPU", choices = NULL),
                              icon("question-circle", id = "help_conn_post", class = "help-icon")
                          ),
                          fluidRow(
                            column(4, numericInput("conn_weight", "Initial Weight", value = 0.1, min = 0, max = 1, step = 0.1)),
                            column(4, numericInput("conn_alpha", "Increment Rate (α)", value = 0.5, min = 0, max = 1, step = 0.1)),
                            column(4, numericInput("conn_beta", "Decrement Rate (β)", value = 0.1, min = 0, max = 1, step = 0.1))
                          ),
                          fluidRow(
                            column(6, numericInput("conn_alpha_prime", "Inhibitory Increment Rate (α')", value = 0.5, min = 0, max = 1, step = 0.1)),
                            column(6, numericInput("conn_beta_prime", "Inhibitory Decrement Rate (β')", value = 0.1, min = 0, max = 1, step = 0.1))
                          ),
                          fluidRow(
                            column(6, 
                                   actionButton("add_connection", "Add Connection", style = "width: 100%;"),
                                   icon("question-circle", id = "help_add_connection", class = "help-icon")
                            ),
                            column(6,
                                   div(
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
                          div(style = "display: flex; align-items: center;",
                              textInput("conn_file_name", "Connections File Name"),
                              icon("question-circle", id = "help_conn_file_name", class = "help-icon")
                          ),
                          div(style = "display: flex; align-items: center;",
                              textInput("conn_file_path", "Connections Directory Path"),
                              icon("question-circle", id = "help_conn_file_path", class = "help-icon")
                          ),
                          actionButton("save_connections", "Save Connections"),
                          actionButton("import_connections", "Import Connections")
                        )
                      ),
                      fluidRow(
                        box(
                          title = "Network Visualization",
                          width = 12,
                          visNetworkOutput("network_plot", height = "400px"),
                          actionButton("reorganize_network", "Reorganize Network")
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
                          div(style = "display: flex; align-items: center;",
                              actionButton("create_architecture", "Create Architecture")
                          )
                        )
                      )
              ),
              tabItem(tabName = "trials",
                      fluidRow(
                        box(
                          title = "Trial Configuration",
                          width = 12,
                          div(style = "display: flex; align-items: center;",
                              textInput("trial_name", "Trial Type Name"),
                              icon("question-circle", id = "help_trial_name", class = "help-icon")
                          ),
                          div(style = "display: flex; align-items: center;",
                              numericInput("num_moments", "Number of Time Moments", 
                                           value = 1, min = 1, max = 10, step = 1),
                              icon("question-circle", id = "help_num_moments", class = "help-icon")
                          ),
                          uiOutput("dynamic_inputs"),
                          checkboxInput("learning_rule", "Active Learning Rule", value = TRUE),
                          actionButton("add_trial", "Add Trial"),
                          actionButton("add_iti", "Add ITI")
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
                          div(style = "display: flex; align-items: center;",
                              textInput("trials_file_name", "Trials File Name (without extension)"),
                              icon("question-circle", id = "help_trials_file_name", class = "help-icon")
                          ),
                          div(style = "display: flex; align-items: center;",
                              textInput("trials_file_path", "Trials Directory Path"),
                              icon("question-circle", id = "help_trials_file_path", class = "help-icon")
                          ),
                          actionButton("save_trials", "Save Trials"),
                          actionButton("import_trials", "Import Trials")
                        )
                      )
              ),
              tabItem(tabName = "contingencies",
                      fluidRow(
                        box(
                          title = "Configure Contingencies",
                          width = 12,
                          div(style = "display: flex; align-items: center;",
                              textInput("phase_name", "Phase or Condition Name"),
                              icon("question-circle", id = "help_phase_name", class = "help-icon")
                          ),
                          div(style = "display: flex; align-items: center;",
                              selectInput("presentation_mode", "Trial Presentation Mode",
                                          choices = c("Random", "In bulk", "Alternated")),
                              icon("question-circle", id = "help_presentation_mode", class = "help-icon")
                          ),
                          div(style = "display: flex; align-items: center;",
                              selectInput("trial_types", "Trial Types to Present",
                                          choices = NULL, multiple = TRUE),
                              icon("question-circle", id = "help_trial_types", class = "help-icon")
                          ),
                          uiOutput("trial_numbers"),
                          checkboxInput("reset_activations", "Reset Activations", value = TRUE),
                          conditionalPanel(
                            condition = "!input.reset_activations",
                            numericInput("min_iti", "Minimum ITI Value", value = 30, min = 1),
                            numericInput("max_iti", "Maximum ITI Value", value = 30, min = 1),
                            div(style = "display: flex; align-items: center;",
                                selectInput("iti_trial", "Add Created ITI", choices = NULL),
                                icon("question-circle", id = "help_iti_trial", class = "help-icon")
                            )
                          ),
                          actionButton("add_contingency", "Add Contingency")
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
                          div(style = "display: flex; align-items: center;",
                              textInput("contingencies_file_name", "Contingencies File Name"),
                              icon("question-circle", id = "help_contingencies_file_name", class = "help-icon")
                          ),
                          div(style = "display: flex; align-items: center;",
                              textInput("contingencies_file_path", "Contingencies Directory Path"),
                              icon("question-circle", id = "help_contingencies_file_path", class = "help-icon")
                          ),
                          actionButton("save_contingencies", "Save Contingencies"),
                          actionButton("import_contingencies", "Import Contingencies")
                        )
                      )
              ),
              tabItem(tabName = "simulate",
                      fluidRow(
                        box(
                          title = "Run Simulation",
                          width = 12,
                          div(style = "display: flex; align-items: center;",
                              numericInput("num_simulations", "Number of Networks", value = 1, min = 1, step = 1),
                              icon("question-circle", id = "help_num_simulations", class = "help-icon")
                          ),
                          actionButton("run_simulation", "Run Simulation"),
                          actionButton("save_networks", "Save Networks"),
                          div(style = "display: flex; align-items: center;",
                              textInput("sim_file_name", "Simulation File Name"),
                              icon("question-circle", id = "help_sim_file_name", class = "help-icon")
                          ),
                          div(style = "display: flex; align-items: center;",
                              textInput("sim_file_path", "Simulation Directory Path"),
                              icon("question-circle", id = "help_sim_file_path", class = "help-icon")
                          ),
                          actionButton("save_simulation", "Save Simulation"),
                          actionButton("load_simulation", "Load Simulation"),
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
              tabItem(tabName = "results_individual",
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
                          actionButton("graficar", "Plot Results")
                        )
                      ),
                      fluidRow(
                        box(
                          title = "Activation Plot",
                          width = 12,
                          selectizeInput("activations_units", "Select Units", 
                                         choices = NULL, multiple = TRUE),
                          selectInput("selected_timestep", "Select Time Step", choices = NULL),
                          plotlyOutput("activations_plot", height = "400px")
                        )
                      ),
                      fluidRow(
                        box(
                          title = "Connection Weights Plot",
                          width = 12,
                          selectizeInput("weights_connections", "Select Connections", 
                                         choices = NULL, multiple = TRUE),
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
                                      choices = c("Mean" = "mean", "Median" = "median")),
                          selectInput("aggregate_error", "Select Error", 
                                      choices = c("Standard Error" = "se", "Standard Deviation" = "sd")),
                          plotlyOutput("aggregate_plot", height = "400px")
                        )
                      )
              ),
              tabItem(tabName = "results_general",
                      fluidRow(
                        box(
                          title = "Plot General Results",
                          width = 12,
                          actionButton("graficar_general", "Plot General Results")
                        )
                      ),
                      fluidRow(
                        box(
                          title = "Activation Plot (All Networks)",
                          width = 12,
                          selectizeInput("activations_units_general", "Select Units", 
                                         choices = NULL, multiple = TRUE),
                          selectInput("selected_timestep_general", "Select Time Step", choices = NULL),
                          plotlyOutput("activations_plot_general", height = "600px")
                        )
                      ),
                      fluidRow(
                        box(
                          title = "Connection Weights Plot (All Networks)",
                          width = 12,
                          selectizeInput("weights_connections_general", "Select Connections", 
                                         choices = NULL, multiple = TRUE),
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
                                      choices = c("Mean" = "mean", "Median" = "median")),
                          selectInput("aggregate_error_general", "Select Error", 
                                      choices = c("Standard Error" = "se", "Standard Deviation" = "sd")),
                          plotlyOutput("aggregate_plot_general", height = "400px")
                        )
                      )
              ),
              tabItem(tabName = "help",
                      fluidRow(
                        box(
                          title = "User Guide",
                          width = 12,
                          tabsetPanel(
                            tabPanel("Introduction",
                                     h4("Welcome to the DDM Simulator"),
                                     p("This guide will help you effectively use the Diffuse Discrepancy Model (DDM) Simulator."),
                                     h4("Model features:"),
                                     tags$ul(
                                       tags$li("Based on principles of neuroanatomy and neurophysiology."),
                                       tags$li("Uses activation and learning rules to model the behavior of neural processing units (NPUs)."),
                                       tags$li("Incorporates hippocampal and dopaminergic systems in the learning process."),
                                       tags$li("Allows simulation of both Pavlovian and operant conditioning.")
                                     )
                            ),
                            tabPanel("Network Architecture",
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
                                     div(style = "overflow-x: auto;",
                                         tags$table(class = "table table-bordered table-striped table-hover",
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
                                     div(style = "overflow-x: auto;",
                                         tags$table(class = "table table-bordered table-striped table-hover",
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
                                                        tags$td("0.5"), tags$td("0.1"), tags$td("0.5"), tags$td("0.1")
                                                      ),
                                                      tags$tr(
                                                        tags$td("S.1"), tags$td("S.2"), tags$td("0.1"),
                                                        tags$td("0.5"), tags$td("0.1"), tags$td("0.5"), tags$td("0.1")
                                                      ),
                                                      tags$tr(
                                                        tags$td("S.2"), tags$td("H1"), tags$td("0.1"),
                                                        tags$td("0.5"), tags$td("0.1"), tags$td("0.5"), tags$td("0.1")
                                                      )
                                                    )
                                         )
                                     ),
                                     p("Note: The NPU names are examples. You can use any name you want for your units.")
                            ),
                            tabPanel("Create Trials",
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
                            tabPanel("Contingency Configuration",
                                     h4("Contingency Design"),
                                     p("In this section, you can configure the experimental contingencies:"),
                                     tags$ul(
                                       tags$li("Define experimental phases."),
                                       tags$li("Specify the trial presentation mode (Random, In block, Alternated)."),
                                       tags$li("Configure inter-trial intervals (ITI) if necessary.")
                                     ),
                                     p("The contingencies determine how different types of trials are presented during the simulation."),
                                     p("Typically, Random is used for training trials and In block for test trials."),
                                     h4("Example of Contingencies file (CSV):"),
                                     div(style = "overflow-x: auto;",
                                         tags$table(class = "table table-bordered table-striped table-hover",
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
                            tabPanel("Simulation",
                                     h4("Running Simulations"),
                                     p("In this section, you can run your simulations:"),
                                     tags$ul(
                                       tags$li("Specify the number of networks to simulate."),
                                       tags$li("Start the simulation and monitor its progress."),
                                       tags$li("Save and load simulations for later analysis.")
                                     ),
                                     p("During the simulation, the model updates the NPU activations and connection weights according to the activation and learning rules.")
                            ),
                            tabPanel("Results Analysis",
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
                          p("For more information, visit the experimental and theoretical research laboratory in Learning, Conditioning and Adaptive Behavior: ",
                            a("Website", href = "http://www.ceic.cucba.udg.mx/Investigacion/laboratorios?id=13", target = "_blank"))
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
    # Tooltips
    bsTooltip("help_npe_name", "Name the network unit. It is recommended NOT to use apostrophes or quotes, and to use a short name. For example: S.1", placement = "right", trigger = "hover"),
    bsTooltip("help_npe_type", "Select whether the unit is excitatory or inhibitory", placement = "right", trigger = "hover"),
    bsTooltip("help_npe_layer", "Select the layer to which the unit belongs. The dopaminergic unit and US already exist by default", placement = "right", trigger = "hover"),
    bsTooltip("help_npe_file_name", "Name of the file where the NPUs will be saved", placement = "right", trigger = "hover"),
    bsTooltip("help_npe_file_path", "Path of the directory where the NPUs file will be saved. Manually copy the path to save the file", placement = "right", trigger = "hover"),
    bsTooltip("help_conn_pre", "Select the source NPU of the connection", placement = "right", trigger = "hover"),
    bsTooltip("help_conn_post", "Select the target NPU of the connection", placement = "right", trigger = "hover"),
    bsTooltip("help_add_connection", "Add the defined connection to the network", placement = "right", trigger = "hover"),
    bsTooltip("help_conn_file_name", "Name of the file where the connections will be saved", placement = "right", trigger = "hover"),
    bsTooltip("help_conn_file_path", "Path of the directory where the connections file will be saved. Manually copy the path to save the file", placement = "right", trigger = "hover"),
    bsTooltip("help_trial_name", "Identifying name for the trial type", placement = "right", trigger = "hover"),
    bsTooltip("help_num_moments", "Number of time moments in the trial", placement = "right", trigger = "hover"),
    bsTooltip("help_trials_file_name", "Name of the file where the trials will be saved", placement = "right", trigger = "hover"),
    bsTooltip("help_trials_file_path", "Path of the directory where the trials file will be saved", placement = "right", trigger = "hover"),
    bsTooltip("help_phase_name", "Name of the experimental phase or condition", placement = "right", trigger = "hover"),
    bsTooltip("help_presentation_mode", "Mode of presentation of the trials", placement = "right", trigger = "hover"),
    bsTooltip("help_trial_types", "Select the types of trials to present in this phase", placement = "right", trigger = "hover"),
    bsTooltip("help_iti_trial", "Select the ITI trial to add", placement = "right", trigger = "hover"),
    bsTooltip("help_contingencies_file_name", "Name of the file where the contingencies will be saved", placement = "right", trigger = "hover"),
    bsTooltip("help_contingencies_file_path", "Path of the directory where the contingencies file will be saved", placement = "right", trigger = "hover"),
    bsTooltip("help_num_simulations", "Number of networks to simulate", placement = "right", trigger = "hover"),
    bsTooltip("help_sim_file_name", "Name of the file where the simulation will be saved", placement = "right", trigger = "hover"),
    bsTooltip("help_sim_file_path", "Path of the directory where the simulation file will be saved", placement = "right", trigger = "hover")
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
      x = runif(2),
      y = runif(2),
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
    
    # Añadir NPE
    observeEvent(input$add_npe, {
      new_npe <- data.frame(
        NPE = input$npe_name,
        Type = input$npe_type,
        Layer = input$npe_layer,
        Activation = input$npe_activation,
        Temporal.Summation = input$npe_temporal_summation,
        Activation.Decay = input$npe_decay,
        mu = input$npe_mu,
        sigma = input$npe_sigma,
        logisSigma = input$npe_logistic_slope,
        x = runif(1),
        y = runif(1),
        stringsAsFactors = FALSE
      )
      current_npes <- rbind(npes(), new_npe)
      npes(current_npes)
      updateSelectInput(session, "conn_pre", choices = current_npes$NPE)
      updateSelectInput(session, "conn_post", choices = current_npes$NPE)
    })
    
    # Añadir Conexión
    observeEvent(input$add_connection, {
      req(input$conn_pre != input$conn_post)
      new_connection <- data.frame(
        PreSinapticNPE = input$conn_pre,
        PostSinapticNPE = input$conn_post,
        Weight = input$conn_weight,
        alpha = input$conn_alpha,
        beta = input$conn_beta,
        alpha_prime = input$conn_alpha_prime,
        beta_prime = input$conn_beta_prime,
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
      
      nodes <- npes() %>%
        mutate(
          id = NPE,
          label = NPE,
          shape = case_when(
            Layer == "PrimarySensory" ~ "square",
            Layer == "US" ~ "hexagon",
            Type == "Inhibitory" ~ "triangle",
            TRUE ~ "circle"
          ),
          color.border = case_when(
            Layer == "PrimarySensory" ~ "black",
            Layer == "AssociativeSensory" ~ "orange",
            Layer == "Hippocampal" ~ "green",
            Layer == "AssociativeMotor" ~ "lightblue",
            Layer == "PrimaryMotor" ~ "blue",
            Layer == "US" ~ "red",
            Layer == "Dopaminergic" ~ "purple",
            TRUE ~ "gray"
          ),
          size = case_when(
            Layer == "PrimarySensory" ~ 16,
            Layer == "US" ~ 16,
            TRUE ~ 25
          ),
          font = list(size = 16),
          title = paste("NPE:", NPE, "<br>Type:", Type, "<br>Layer:", Layer)  # Tooltip
        )
      
      edges <- connections() %>%
        mutate(
          from = PreSinapticNPE,
          to = PostSinapticNPE,
          arrows = "to",
          color = case_when(
            Weight == 1 ~ "black",
            PreSinapticNPE %in% (nodes %>% filter(Type == "Inhibitory") %>% pull(NPE)) ~ "red",
            TRUE ~ "gray"
          ),
          width = Weight * 2,
          title = paste("Weight:", Weight, "<br>alpha:", alpha, "<br>beta:", beta)  # Tooltip
        )
      
      visNetwork(nodes, edges, width = "100%", height = "400px") %>%
        visOptions(highlightNearest = list(enabled = TRUE, degree = 1)) %>%
        visPhysics(stabilization = FALSE) %>%
        visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) %>%
        visLayout(randomSeed = 123) %>%
        visEvents(type = "once", startStabilizing = "function() {
          this.moveTo({scale:0.8})
        }")
    })
    
    # Actualizar posiciones de los nodos
    observe({
      network_data <- input$network_plot_positions
      if (!is.null(network_data)) {
        node_positions <- network_data$nodes %>% 
          select(id, x, y)
        updated_npes <- npes() %>% 
          left_join(node_positions, by = c("NPE" = "id")) %>% 
          mutate(
            x = coalesce(x, 0),
            y = coalesce(y, 0)
          )
        npes(updated_npes)
      }
    })
    
    # Reorganizar la red
    observeEvent(input$reorganize_network, {
      npes_updated <- npes() %>%
        mutate(
          x = runif(n()),
          y = runif(n())
        )
      npes(npes_updated)
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
      req(input$npe_file_name != "", input$npe_file_path != "")
      tryCatch({
        filename <- ensure_csv_extension(input$npe_file_name)
        full_path <- file.path(input$npe_file_path, filename)
        write.csv(npes(), full_path, row.names = FALSE)
        showNotification("NPEs saved successfully", type = "message")
      }, error = function(e) {
        showNotification(paste("Error saving NPEs:", e$message), type = "error", duration = NULL)
      })
    })
    
    # Guardar Conexiones
    observeEvent(input$save_connections, {
      req(input$conn_file_name != "", input$conn_file_path != "")
      tryCatch({
        filename <- ensure_csv_extension(input$conn_file_name)
        full_path <- file.path(input$conn_file_path, filename)
        write.csv(connections(), full_path, row.names = FALSE)
        showNotification("Connections saved successfully", type = "message")
      }, error = function(e) {
        showNotification(paste("Error saving Connections:", e$message), type = "error", duration = NULL)
      })
    })
    
    # Importar NPEs
    observeEvent(input$import_npes, {
      req(input$npe_file_name != "", input$npe_file_path != "")
      tryCatch({
        filename <- ensure_csv_extension(input$npe_file_name)
        full_path <- file.path(input$npe_file_path, filename)
        if (!file.exists(full_path)) {
          stop(paste("The file does not exist:", full_path))
        }
        imported_npes <- as.data.frame(read_csv(full_path))
        npes(imported_npes)
        updateSelectInput(session, "conn_pre", choices = imported_npes$NPE)
        updateSelectInput(session, "conn_post", choices = imported_npes$NPE)
        showNotification("NPEs imported successfully", type = "message")
      }, error = function(e) {
          showNotification(paste("Error importing NPEs:", e$message), type = "error", duration = NULL)
      })
    })
    
    # Importar Conexiones
    observeEvent(input$import_connections, {
      req(input$conn_file_name != "", input$conn_file_path != "")
      tryCatch({
        filename <- ensure_csv_extension(input$conn_file_name)
        full_path <- file.path(input$conn_file_path, filename)
        if (!file.exists(full_path)) {
          stop(paste("The file does not exist:", full_path))
        }
        imported_connections <- as.data.frame(read_csv(full_path))
        connections(imported_connections)
        showNotification("Connections imported successfully", type = "message")
      }, error = function(e) {
        showNotification(paste("Error importing Connections:", e$message), type = "error", duration = NULL)
      })
    })
    
    # Observador para el número de momentos temporales
    observeEvent(input$num_moments, {
      if (is.null(input$num_moments) || is.na(input$num_moments) || input$num_moments < 1) {
        updateNumericInput(session, "num_moments", value = 1)
        showNotification("The minimum number of time steps is 1", type = "warning")
      }
    }, ignoreInit = TRUE, ignoreNULL = FALSE)
    
    # Observadores para limitar los valores de entrada a 0-1
    observe({
      req(input$num_moments)
      lapply(1:input$num_moments, function(moment) {
        lapply(npes()$NPE[npes()$Layer %in% c("US", "PrimarySensory")], function(npe) {
          input_id <- paste0("input_", npe, "_", moment)
          observeEvent(input[[input_id]], {
            if (is.null(input[[input_id]]) || is.na(input[[input_id]]) || input[[input_id]] < 0) {
              updateNumericInput(session, input_id, value = 0)
              showNotification("Values must be between 0 and 1", type = "warning")
            } else if (input[[input_id]] > 1) {
              updateNumericInput(session, input_id, value = 1)
              showNotification("Values must be between 0 and 1", type = "warning")
            }
          }, ignoreInit = TRUE, ignoreNULL = FALSE)
        })
      })
    })
    
    # Regenerar inputs dinámicos cuando cambia el número de momentos
    observeEvent(input$num_moments, {
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
    }, ignoreInit = TRUE)
    
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
          modalButton("Cancelar"),
          actionButton("save_iti", "Save ITI")
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
                                          step = 0.1)
      }
      input_list$iti_learning_rule <- checkboxInput("iti_learning_rule", "Regla de aprendizaje activa", value = TRUE)
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
                        choices = names(current_trials)[sapply(current_trials, function(x) is.character(x) && length(x) == 1)])
    })
    
    # Mostrar botones de ensayos
    output$trial_buttons <- renderUI({
      all_trials <- trials()
      if (length(all_trials) == 0) return(NULL)
      
      tagList(
        lapply(names(all_trials), function(trial_type) {
          div(
            actionButton(inputId = paste0("trial_", trial_type), 
                         label = trial_type, 
                         class = "btn-block"),
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
          lapply(seq(1, length(trial_parts)-1, 2), function(i) {
            p(paste(trial_parts[i], ":", trial_parts[i+1]))
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
              lapply(seq(1, length(trial_parts)-1, 2), function(j) {
                p(paste(trial_parts[j], ":", trial_parts[j+1]))
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
                     value = 1, min = 1)
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
      if (length(cont) == 0) return(NULL)
      
      cont_df <- data.frame(
        Contingencia = cont,
        stringsAsFactors = FALSE
      )
      
      datatable(cont_df, options = list(pageLength = 10, scrollX = TRUE, scrollY = "300px"))
    })
    
    # Observador para el número de simulaciones
    observeEvent(input$num_simulations, {
      if (is.na(input$num_simulations) || input$num_simulations < 1) {
        updateNumericInput(session, "num_simulations", value = 1)
        showNotification("The minimum number of simulations is 1", type = "warning")
      }
    }, ignoreInit = TRUE)
    
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
        showNotification(paste("Missing necessary data for the simulation:", 
                               paste(missing_components, collapse = ", "), 
                               ". Ensure that you have created the architecture and defined the trials and contingencies."), 
                         type = "error", duration = NULL)
        return()
      }
      
      output$simulation_status <- renderText("Simulation in progress...")
      updateProgressBar(session, "sim_progress", value = 0)
      
      all_results <- list()
      
      for(i in 1:input$num_simulations) {
        tryCatch({
          # Usar los datos del entorno global
          NPEs <- get("NPEs", envir = .GlobalEnv)
          Connections <- get("Connections", envir = .GlobalEnv)
          contingencies <- get("contingencies", envir = .GlobalEnv)
          trials <- get("trials", envir = .GlobalEnv)
          
          # Calcular HasITI
          HasITI <- sapply(contingencies, function(x) as.logical(trimws(unlist(strsplit(x,",")))[5]))
          assign("HasITI", HasITI, envir = .GlobalEnv)
          
          # Crear TimeSteps
          TimeSteps <- Create.Phases(contingencies, trials)
          assign("TimeSteps", TimeSteps, envir = .GlobalEnv)
          
          # Ejecutar la simulación
          datos <- Simulate.DBP(NPEs, Connections, TimeSteps, HasITI)
          
          # Guardar los resultados
          all_results[[i]] <- datos
          
          updateProgressBar(session, "sim_progress", value = (i / input$num_simulations) * 100)
          output$simulation_status <- renderText(paste("Simulation", i, "of", input$num_simulations, "completed"))
        }, error = function(e) {
          output$simulation_status <- renderText(paste("Error in simulation", i, ":", e$message))
        })
      }
      
      # Guardar todos los resultados
      simulation_results(all_results)
      assign("datos_simulacion", all_results, envir = .GlobalEnv)
      output$simulation_status <- renderText("All simulations completed")
      
      # Actualizar el selector de simulaciones
      updateSelectInput(session, "selected_simulation", 
                        choices = paste("Network", 1:length(all_results)))
    })
    
    # Guardar simulación
    observeEvent(input$save_simulation, {
      req(input$sim_file_name != "", input$sim_file_path != "")
      
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
        showNotification(paste("The simulation cannot be saved. The following components are missing:", 
                               paste(missing_components, collapse = ", ")), 
                         type = "error", duration = NULL)
        return()
      }
      
      simulation_data <- list(
        NPEs = get("NPEs", envir = .GlobalEnv),
        Connections = get("Connections", envir = .GlobalEnv),
        trials = get("trials", envir = .GlobalEnv),
        contingencies = get("contingencies", envir = .GlobalEnv),
        datos_simulacion = get("datos_simulacion", envir = .GlobalEnv)
      )
      
      tryCatch({
        filename <- ensure_rds_extension(input$sim_file_name)
        full_path <- file.path(input$sim_file_path, filename)
        saveRDS(simulation_data, file = full_path)
        showNotification("Simulation saved successfully", type = "message")
      }, error = function(e) {
        showNotification(paste("Error saving the simulation:", e$message), type = "error", duration = NULL)
      })
    })
    
    # Cargar simulación
    observeEvent(input$load_simulation, {
      req(input$sim_file_name != "", input$sim_file_path != "")
      tryCatch({
        filename <- ensure_rds_extension(input$sim_file_name)
        full_path <- file.path(input$sim_file_path, filename)
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
                          choices = names(simulation_data$trials)[sapply(simulation_data$trials, function(x) is.character(x) && length(x) == 1)])
        
        # Actualizar el selector de simulaciones
        updateSelectInput(session, "selected_simulation", 
                          choices = paste("Network", 1:length(simulation_data$datos_simulacion)))
        
        showNotification("Simulation successfully loaded", type = "message")
      }, error = function(e) {
        showNotification(paste("Error loading the simulation:", e$message), type = "error", duration = NULL)
      })
    })
    
    # Mostrar resultados de la simulación
    output$simulation_results <- renderDT({
      req(simulation_results())
      selected_sim <- as.numeric(gsub("Network ", "", input$selected_simulation))
      datatable(head(simulation_results()[[selected_sim]], 1000), 
                options = list(scrollX = TRUE, scrollY = "400px"))
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
          npes_in_trial <- unique(configuration_step[seq(1, length(configuration_step)-1, 2)])
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
                    choices = c("All networks in one file", "Individual networks")),
        textInput("save_networks_filename", "File name (without extension):"),
        textInput("save_networks_path", "Save path:"),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_save_networks", "Save")
        )
      ))
    })
    
    observeEvent(input$confirm_save_networks, {
      req(input$save_networks_filename, input$save_networks_path)
      tryCatch({
        results <- simulation_results()
        if (input$save_option == "All networks in one file") {
          combined_results <- do.call(rbind, lapply(seq_along(results), function(i) {
            cbind(Network = i, results[[i]])
          }))
          filename <- file.path(input$save_networks_path, paste0(input$save_networks_filename, ".csv"))
          write.csv(combined_results, filename, row.names = FALSE)
        } else {
          for (i in seq_along(results)) {
            filename <- file.path(input$save_networks_path, paste0(input$save_networks_filename, "_Network_", i, ".csv"))
            write.csv(results[[i]], filename, row.names = FALSE)
          }
        }
        removeModal()
        showNotification("Networks saved successfully", type = "message")
      }, error = function(e) {
        showNotification(paste("Error saving the networks:", e$message), type = "error")
      })
    })
    
    # Procesamiento de datos para gráficos
    processed_data <- reactiveVal(NULL)
    
    observeEvent(input$graficar, {
      if (!exists("datos_simulacion", envir = .GlobalEnv) || length(get("datos_simulacion", envir = .GlobalEnv)) == 0) {
        showNotification("There is no simulation data to plot. Please run a simulation first.", type = "error")
        return()
      }
      
      all_data <- get("datos_simulacion", envir = .GlobalEnv)
      selected_sim <- as.numeric(gsub("Network ", "", input$selected_simulation))
      data <- all_data[[selected_sim]]
      
      # Asegurarse de que 'Phase' sea character y 'Trial' y 'TimeStep' sean numeric
      data$Phase <- as.character(data$Phase)
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
                           selected = ordered_units[npes_data$Layer[match(ordered_units, npes_data$NPE)] == "PrimaryMotor"])
      
      # Identificar unidades PrimaryMotor
      primary_motor_units <- ordered_units[npes_data$Layer[match(ordered_units, npes_data$NPE)] == "PrimaryMotor"]
      
      updateSelectizeInput(session, "aggregate_units",
                           choices = ordered_units,
                           selected = primary_motor_units)
      updateSelectInput(session, "aggregate_phase",
                        choices = unique(activations$Phase))
      updateSelectInput(session, "selected_timestep",
                        choices = unique(activations$TimeStep),
                        selected = max(unique(activations$TimeStep)) - 1)
      updateSelectInput(session, "aggregate_timestep",
                        choices = unique(activations$TimeStep),
                        selected = max(unique(activations$TimeStep)) - 1)
      
      # Actualizar opciones de conexiones
      updateSelectizeInput(session, "weights_connections", 
                           choices = unique(weights$Connection),
                           selected = unique(weights$Connection)[1:min(5, length(unique(weights$Connection)))])
    })
    
    # Gráfico de Activaciones (Resultados individuales)
    output$activations_plot <- renderPlotly({
      req(processed_data(), input$activations_units, input$selected_timestep)
      data <- processed_data()$activations %>%
        filter(Unit %in% input$activations_units, TimeStep == input$selected_timestep)
      
      p <- ggplot(data, aes(x = Trial, y = Activation, color = Unit)) +
        geom_line(size = 0.6) +  # Línea más delgada
        facet_wrap(~Phase, scales = "free_x") +
        theme_minimal() +
        labs(title = paste("Activations per trial (Timestep", input$selected_timestep, ")"),
             x = "Trial", y = "Activation")
      
      ggplotly(p) %>%
        layout(legend = list(orientation = "h", y = -0.2)) %>%
        config(displayModeBar = TRUE) %>%
        add_trace(
          text = ~paste(
            "Trial:", Trial,
            "<br>Activación:", round(Activation, 5),
            "<br>Unidad:", Unit
          ),
          hoverinfo = "text"
        )
    })
    
    # Gráfico de Pesos de Conexiones (Resultados individuales)
    output$weights_plot <- renderPlotly({
      req(processed_data(), input$weights_connections, input$selected_timestep)
      data <- processed_data()$weights %>%
        filter(Connection %in% input$weights_connections, TimeStep == input$selected_timestep)
      
      if(nrow(data) == 0) {
        return(plot_ly() %>% add_annotations(text = "There is no data to display", showarrow = FALSE))
      }
      
      p <- ggplot(data, aes(x = Trial, y = Weight, color = Connection)) +
        geom_line(size = 0.8) +  # Línea más delgada
        facet_wrap(~Phase, scales = "free_x") +
        theme_minimal() +
        labs(title = paste("Connection weights per trial (Timestep", input$selected_timestep, ")"),
             x = "Trial", y = "Weights")
      
      ggplotly(p) %>%
        layout(legend = list(orientation = "h", y = -0.2)) %>%
        config(displayModeBar = TRUE) %>%
        add_trace(
          text = ~paste(
            "Trial:", Trial,
            "<br>Weight:", round(Weight, 5),
            "<br>Connection:", Connection
          ),
          hoverinfo = "text"
        )
    })
    
    # Gráfico de Medidas Agregadas
    output$aggregate_plot <- renderPlotly({
      req(processed_data(), input$aggregate_phase, input$aggregate_units, 
          input$aggregate_timestep, input$aggregate_measure, input$aggregate_error)
      
      data <- processed_data()$activations %>%
        filter(Phase == input$aggregate_phase,
               TimeStep == as.numeric(input$aggregate_timestep),
               Unit %in% input$aggregate_units)
      
      if(nrow(data) == 0) {
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
        add_trace(data = aggregated_data, x = ~Unit, y = as.formula(paste0("~", y_value)), type = "bar",
                  marker = list(color = 'rgba(158,202,225,0.6)', line = list(color = 'rgb(8,48,107)', width = 1.5)),
                  error_y = list(type = "data", array = aggregated_data[[error_value]], visible = TRUE),
                  name = "Average") %>%
        layout(title = paste("Measure used:", input$aggregate_measure),
               xaxis = list(title = "Unit"),
               yaxis = list(title = "Average of activations", range = c(0, 1)),
               showlegend = TRUE)
      
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
      
      # Asegurarse de que 'Phase' sea character y 'Trial' y 'TimeStep' sean numeric
      combined_data$Phase <- as.character(combined_data$Phase)
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
                           selected = ordered_units[npes_data$Layer[match(ordered_units, npes_data$NPE)] == "PrimaryMotor"])
      
      # Identificar unidades PrimaryMotor
      primary_motor_units <- ordered_units[npes_data$Layer[match(ordered_units, npes_data$NPE)] == "PrimaryMotor"]
      
      updateSelectizeInput(session, "aggregate_units_general",
                           choices = ordered_units,
                           selected = primary_motor_units)
      updateSelectInput(session, "aggregate_phase_general",
                        choices = unique(activations$Phase))
      updateSelectInput(session, "selected_timestep_general",
                        choices = unique(activations$TimeStep),
                        selected = max(unique(activations$TimeStep)) - 1)
      updateSelectInput(session, "aggregate_timestep_general",
                        choices = unique(activations$TimeStep),
                        selected = max(unique(activations$TimeStep)) - 1)
      
      # Actualizar opciones de conexiones
      updateSelectizeInput(session, "weights_connections_general", 
                           choices = unique(weights$Connection),
                           selected = unique(weights$Connection)[1:min(5, length(unique(weights$Connection)))])
    })
    
    # Gráfico de Activaciones (Todas las Redes)
    output$activations_plot_general <- renderPlotly({
      req(processed_data_general(), input$activations_units_general, input$selected_timestep_general)
      data <- processed_data_general()$activations %>%
        filter(Unit %in% input$activations_units_general, TimeStep == input$selected_timestep_general)
      
      # Calcular el promedio de activación por ensayo y unidad
      avg_data <- data %>%
        group_by(Phase, Trial, Unit) %>%
        summarise(Avg_Activation = mean(Activation, na.rm = TRUE),
                  SD_Activation = sd(Activation, na.rm = TRUE),
                  .groups = 'drop')
      
      p <- ggplot(avg_data, aes(x = Trial, y = Avg_Activation, color = Unit)) +
        geom_line(size = 1.2) +
        facet_wrap(~Phase, scales = "free_x") +
        theme_minimal() +
        labs(title = paste("Average activations per trial (Timestep", input$selected_timestep_general, ")"), 
             x = "Trial", y = "Average activation")
      
      ggplotly(p) %>%
        layout(legend = list(orientation = "h", y = -0.2)) %>%
        config(displayModeBar = TRUE) %>%
        # Añadir información detallada para la vista dinámica
        add_trace(
          text = ~paste(
            "Trial:", Trial,
            "<br>Average Activation:", round(Avg_Activation, 5),
            "<br>NPE:", Unit,
            "<br>Standard deviation:", round(SD_Activation, 5)
          ),
          hoverinfo = "text"
        )
    })
    
    # Gráfico de Pesos de Conexiones (Todas las Redes)
    output$weights_plot_general <- renderPlotly({
      req(processed_data_general(), input$weights_connections_general, input$selected_timestep_general)
      data <- processed_data_general()$weights %>%
        filter(Connection %in% input$weights_connections_general, TimeStep == input$selected_timestep_general)
      
      if(nrow(data) == 0) {
        return(plot_ly() %>% add_annotations(text = "No data to display", showarrow = FALSE))
      }
      
      # Calcular el promedio de peso por ensayo y conexión
      avg_data <- data %>%
        group_by(Phase, Trial, Connection) %>%
        summarise(Avg_Weight = mean(Weight, na.rm = TRUE),
                  SD_Weight = sd(Weight, na.rm = TRUE),
                  .groups = 'drop')
      
      p <- ggplot(avg_data, aes(x = Trial, y = Avg_Weight, color = Connection)) +
        geom_line(size = 1.2) +
        facet_wrap(~Phase, scales = "free_x") +
        theme_minimal() +
        labs(title = paste("Average connection weights per trial (Timestep", input$selected_timestep_general, ")"), 
             x = "Trial", y = "Average weight")
      
      ggplotly(p) %>%
        layout(legend = list(orientation = "h", y = -0.2)) %>%
        config(displayModeBar = TRUE) %>%
        # Añadir información detallada para la vista dinámica
        add_trace(
          text = ~paste(
            "Trial:", Trial,
            "<br>Average Weight:", round(Avg_Weight, 5),
            "<br>Connection", Connection,
            "<br>Standard deviation:", round(SD_Weight, 5)
          ),
          hoverinfo = "text"
        )
    })
    
    # Gráfico de Medidas Agregadas (Todas las Redes)
    output$aggregate_plot_general <- renderPlotly({
      req(processed_data_general(), input$aggregate_phase_general, input$aggregate_units_general, 
          input$aggregate_timestep_general, input$aggregate_measure_general, input$aggregate_error_general)
      
      data <- processed_data_general()$activations %>%
        filter(Phase == input$aggregate_phase_general,
               TimeStep == as.numeric(input$aggregate_timestep_general),
               Unit %in% input$aggregate_units_general)
      
      if(nrow(data) == 0) {
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
                         median(Activation, na.rm = TRUE)),
          .groups = "drop"
        )
      
      y_value <- ifelse(input$aggregate_measure_general == "mean", "Mean", "Median")
      error_value <- ifelse(input$aggregate_error_general == "se", "SE", "SD")
      
      # Crear un vector de unidades únicas ordenadas
      unique_units <- unique(aggregated_data$Unit)
      
      p <- plot_ly() %>%
        add_trace(data = aggregated_data, x = ~factor(Unit, levels = unique_units), y = as.formula(paste0("~", y_value)), type = "bar",
                  marker = list(color = 'rgba(158,202,225,0.6)', line = list(color = 'rgb(8,48,107)', width = 1.5)),
                  error_y = list(type = "data", array = aggregated_data[[error_value]], visible = TRUE),
                  name = "Average") %>%
        add_trace(data = individual_data, x = ~factor(Unit, levels = unique_units), y = ~Value, type = "scatter", mode = "markers",
                  marker = list(color = 'white', size = 8, line = list(color = 'black', width = 1)),
                  name = "Individual networks") %>%
        layout(title = paste("Measure used:", input$aggregate_measure_general),
               xaxis = list(title = "Type of NPE", 
                            type = 'category',
                            categoryorder = "array",
                            categoryarray = unique_units),
               yaxis = list(title = "Average of activations", range = c(0, 1)),
               showlegend = TRUE,
               legend = list(orientation = "h", y = -0.2),
               barmode = 'overlay')
      
      return(p)
    })
    
    # Guardar Ensayos
    observeEvent(input$save_trials, {
      req(input$trials_file_name != "", input$trials_file_path != "")
      tryCatch({
        filename <- ensure_rds_extension(input$trials_file_name)
        full_path <- file.path(input$trials_file_path, filename)
        
        # Guardar la lista de ensayos directamente como un objeto RDS
        saveRDS(trials(), file = full_path)
        showNotification("Trials saved successfully", type = "message")
      }, error = function(e) {
        showNotification(paste("Error saving the trials:", e$message), type = "error", duration = NULL)
      })
    })
    
    # Importar Ensayos
    observeEvent(input$import_trials, {
      req(input$trials_file_name != "", input$trials_file_path != "")
      tryCatch({
        filename <- ensure_rds_extension(input$trials_file_name)
        full_path <- file.path(input$trials_file_path, filename)
        if (!file.exists(full_path)) {
          stop(paste("El archivo no existe:", full_path))
        }
        
        # Leer el archivo RDS
        imported_trials <- readRDS(full_path)
        
        trials(imported_trials)
        assign("trials", imported_trials, envir = .GlobalEnv)
        updateSelectInput(session, "trial_types", choices = names(imported_trials))
        showNotification("Trials imported successfully", type = "message")
      }, error = function(e) {
        showNotification(paste("Error importing trials:", e$message), type = "error", duration = NULL)
      })
    })
    
    # Guardar Contingencias
    observeEvent(input$save_contingencies, {
      req(input$contingencies_file_name != "", input$contingencies_file_path != "")
      tryCatch({
        filename <- ensure_csv_extension(input$contingencies_file_name)
        full_path <- file.path(input$contingencies_file_path, filename)
        
        # Convertir el vector de contingencias a un dataframe
        contingencies_df <- data.frame(Contingencia = contingencies(), stringsAsFactors = FALSE)
        
        # Guardar como CSV
        write.csv(contingencies_df, file = full_path, row.names = FALSE)
        showNotification("Contingencies saved successfully", type = "message")
      }, error = function(e) {
        showNotification(paste("Error saving the contingencies:", e$message), type = "error", duration = NULL)
      })
    })
    
    # Importar Contingencias
    observeEvent(input$import_contingencies, {
      req(input$contingencies_file_name != "", input$contingencies_file_path != "")
      tryCatch({
        filename <- ensure_csv_extension(input$contingencies_file_name)
        full_path <- file.path(input$contingencies_file_path, filename)
        if (!file.exists(full_path)) {
          stop(paste("El archivo no existe:", full_path))
        }
        
        # Leer el CSV
        imported_contingencies_df <- read.csv(full_path, stringsAsFactors = FALSE)
        
        # Convertir el dataframe de vuelta a un vector
        imported_contingencies <- imported_contingencies_df$Contingencia
        
        contingencies(imported_contingencies)
        assign("contingencies", imported_contingencies, envir = .GlobalEnv)
        showNotification("Contingencies imported successfully", type = "message")
      }, error = function(e) {
        showNotification(paste("Error importing contingencies:", e$message), type = "error", duration = NULL)
      })
    })
    
    # Botón de cierre
    observeEvent(input$close_app, {
      showModal(modalDialog(
        title = "Confirm closure",
        "Are you sure you want to close the program?",
        footer = tagList(
          modalButton("Cancelar"),
          actionButton("confirm_close", "Yes, close", class = "btn-danger")
        )
      ))
    })
    
    observeEvent(input$confirm_close, {
      stopApp()
    })
    
    # Observadores para los botones de ayuda
    observeEvent(input$help_npe_name, {
      showModal(modalDialog(
        title = "Ayuda: Nombre del NPE",
        "NName the network unit. It is recommended NOT to use apostrophes or quotes, and to use a short name. For example: S.1"
      ))
    })
    
    observeEvent(input$help_npe_type, {
      showModal(modalDialog(
        title = "Ayuda: Tipo de NPE",
        "Define whether your network unit is excitatory or inhibitory."
      ))
    })
    
    observeEvent(input$help_npe_layer, {
      showModal(modalDialog(
        title = "Ayuda: Capa",
        "Type of unit you are going to create. The dopaminergic unit and US are already created by default."
      ))
    })
    
    # Añade este código junto con los otros observadores de ayuda
    observeEvent(input$reminder_us_d, {
      showTooltip("reminder_us_d", "Remember that the connection between US and D must be with a maximum weight of 1.", placement = "right", trigger = "hover")
    })
    
    observeEvent(input$help_npe_file_name, {
      showModal(modalDialog(
        title = "Ayuda: Nombre del archivo NPEs",
        "Register a name for your created units. This file will be saved with this name."
      ))
    })
    
    observeEvent(input$help_npe_file_path, {
      showModal(modalDialog(
        title = "Ayuda: Ruta del directorio NPEs",
        "Add the path where you want your file to be saved. You can copy and paste it."
      ))
    })
    
    observeEvent(input$help_conn_pre, {
      showModal(modalDialog(
        title = "Ayuda: NPE de origen",
        "Select the source NPE to connect."
      ))
    })
    
    observeEvent(input$help_conn_post, {
      showModal(modalDialog(
        title = "Ayuda: NPE de destino",
        "Select the destination NPE to connect."
      ))
    })
    
    observeEvent(input$help_add_connection, {
      showModal(modalDialog(
        title = "Ayuda: Añadir conexión",
        "
Save the connections created."
      ))
    })
    
    observeEvent(input$help_conn_file_name, {
      showModal(modalDialog(
        title = "Ayuda: Nombre del archivo Conexiones",
        "Register a name for your created connections. This file will be saved with this name."
      ))
    })
    
    observeEvent(input$help_conn_file_path, {
      showModal(modalDialog(
        title = "Ayuda: Ruta del directorio Conexiones",
        "Add the path where you want your file to be saved. You can copy and paste it."
      ))
    })
    
    
    observeEvent(input$help_trial_name, {
      showModal(modalDialog(
        title = "Ayuda: Nombre del tipo de ensayo",
        "Record what name you want for your essay."
      ))
    })
    
    observeEvent(input$help_num_moments, {
      showModal(modalDialog(
        title = "Ayuda: Número de momentos temporales",
        "Shows how many time moments each of your S' (primary sensory) units will have."
      ))
    })
    
    observeEvent(input$help_phase_name, {
      showModal(modalDialog(
        title = "Ayuda: Nombre de la fase o condición",
        "Record the name of your phase or condition, for example: Training."
      ))
    })
    
    observeEvent(input$help_presentation_mode, {
      showModal(modalDialog(
        title = "Ayuda: Modo de presentación de ensayos",
        "Muestra qué tipo de presentación tendrá cada ensayo. Normalmente Random se utiliza para ensayos de entrenamiento e In Bulk para ensayos de prueba."
      ))
    })
    
    observeEvent(input$help_trial_types, {
      showModal(modalDialog(
        title = "Ayuda: Tipos de ensayos a presentar",
        "Record which essays will be presented in that contingency. For example: If you created 2 S' named X.1 and X.2 for training, you can select X.1 and X.2, if you created a trial named XY for testing, you can select only XY."
      ))
    })
    
    observeEvent(input$help_iti_trial, {
      showModal(modalDialog(
        title = "Ayuda: Agregar ITI creado",
        "If you added an ITI you can select it here and decide how many time points it will have as a minimum and maximum."
      ))
    })
    
    observeEvent(input$help_num_simulations, {
      showModal(modalDialog(
        title = "Ayuda: Número de redes",
        "Choose the number of networks you want to simulate."
      ))
    })
    
    observeEvent(input$help_sim_file_name, {
      showModal(modalDialog(
        title = "Ayuda: Nombre del archivo de simulación",
        "Choose what name your file will have."
      ))
    })
    
    observeEvent(input$help_sim_file_path, {
      showModal(modalDialog(
        title = "Ayuda: Ruta del directorio de simulación",
        "Choose the path where you want to save your file, this file will have everything you created previously. Essential if you just want to simulate again without repeating the other steps."
      ))
    })
  }

# Ejecutar la aplicación Shiny
shinyApp(ui = ui, server = server)
