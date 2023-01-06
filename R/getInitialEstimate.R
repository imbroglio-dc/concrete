#' Title
#'
#' @param Data data.table
#' @param Model list
#' @param CVFolds : list
#' @param MinNuisance numeric
#' @param TargetEvent numeric vector
#' @param TargetTime numeric vector
#' @param Regime list
#' @param PropScoreBackend character
#' @param HazEstBackend character
#' @param ReturnModels boolean

getInitialEstimate <- function(Data, Model, CVFolds, MinNuisance, TargetEvent, TargetTime, 
                               Regime, PropScoreBackend, HazEstBackend, ReturnModels) {
    Time <- NULL
    TrtVal <- Data[[attr(Data, "Treatment")]]
    TimeVal <- Data[[attr(Data, "EventTime")]]
    Censored <- any(Data[[attr(Data, "EventType")]] <= 0)
    CovDT <- subset(Data, select = attr(Data, "CovNames")[["ColName"]])
    TrtModel <- try(Model[[attr(Data, "Treatment")]])
    options(warn = 0)
    if (inherits(TrtModel, "try-error"))
        stop("TrtModel must currently be specified in the Model argument as a list named ", 
             "as the Treatment variable, ", attr(Data, "Treatment"))
    
    ## Propensity Scores for Regimes of Interest ----
    cat("Trt: ")
    PropScores <- getPropScore(TrtVal = TrtVal, CovDT = CovDT, TrtModel = TrtModel,
                               MinNuisance = MinNuisance, Regime = Regime,
                               PropScoreBackend = PropScoreBackend, CVFolds = CVFolds, TrtLoss = NULL, 
                               ReturnModels = ReturnModels)
    InitFits <- list(attr(PropScores, "TrtFit"))
    names(InitFits)[1] <- attr(Data, "Treatment")
    cat("Done\n")
    
    ## hazards: Events & censoring ----
    ## baseline hazards for obs times + target times ----
    HazTimes <- sort(unique(c(TargetTime, TimeVal)))
    HazTimes <- HazTimes[HazTimes <= max(TargetTime)]
    Hazards <- data.table("Time" = c(0, HazTimes))
    
    cat("Hazards: ")
    HazFits <- getHazFit(Data = Data, Model = Model, CVFolds = CVFolds, Hazards = Hazards, 
                         HazEstBackend = HazEstBackend, ReturnModels = ReturnModels)
    InitFits <- c(InitFits, lapply(HazFits, function(HF) return(attr(HF, "HazSL"))))
    HazSurvPreds <- getHazSurvPred(Data = Data, HazFits = HazFits, MinNuisance = MinNuisance,
                                   TargetEvent = TargetEvent, TargetTime = TargetTime,
                                   Regime = Regime, HazEstBackend)
    InitialEstimates <- lapply(seq_along(PropScores), function(a) {
        if (Censored) {
            NuisanceWeight <- sapply(seq_along(PropScores[[a]]), function(i) {
                PropScores[[a]][i] * HazSurvPreds[[a]][["Survival"]][["LaggedCensSurv"]][, i]
            })   
        } else {
            NuisanceWeight <- matrix(PropScores[[a]], 
                                     nrow = nrow(HazSurvPreds[[a]][["Survival"]][["TotalSurv"]]), 
                                     ncol = ncol(HazSurvPreds[[a]][["Survival"]][["TotalSurv"]]), 
                                     byrow = TRUE) 
        }
        NuisanceWeight <- 1 / truncNuisanceDenom(NuisanceDenom = NuisanceWeight, 
                                                 MinNuisance = MinNuisance, 
                                                 RegimeName = names(PropScores)[a])
        return(list("PropScore" = PropScores[[a]],
                    "Hazards" = HazSurvPreds[[a]][["Hazards"]],
                    "EvntFreeSurv" = HazSurvPreds[[a]][["Survival"]][["TotalSurv"]],
                    "NuisanceWeight" = NuisanceWeight))
    })
    
    names(InitialEstimates) <- names(Regime)
    attr(InitialEstimates, "Times") <- Hazards[["Time"]]
    attr(InitialEstimates, "InitFits") <- InitFits
    cat("Done\n")
    return(InitialEstimates)
}

truncNuisanceDenom <- function(NuisanceDenom, MinNuisance, RegimeName) {
    if (is.function(MinNuisance)) {
        warning("MinNuisance functions are not yet supported.")
        MinNuisance <- 5 / log(ncol(NuisanceDenom)) / sqrt(ncol(NuisanceDenom))
    }
    if (is.numeric(MinNuisance) & length(MinNuisance) == 1) {
        if (MinNuisance < 1 & MinNuisance > 0) {
            if (min(NuisanceDenom) < MinNuisance) {
                PositivityWarning <- paste(
                    "For Intervention \"", RegimeName, "\", ", 
                    round(mean(NuisanceDenom < MinNuisance), 3) * 100, "% of the G-related ", 
                    "nuisance parameter estimates were bounded", 
                    # sum(apply(NuisanceDenom, 2, function(subj) any(subj < MinNuisance))), 
                    # "/", ncol(NuisanceDenom), " subjects", 
                    " to ", signif(MinNuisance, 3), sep = "")
                attr(NuisanceDenom, "original") <- NuisanceDenom
                attr(NuisanceDenom, "message") <- PositivityWarning
                NuisanceDenom[NuisanceDenom < MinNuisance] <- MinNuisance
            } else {
                attr(NuisanceDenom, "message") <- paste(
                    "For Intervention \"", RegimeName, "\", ", 
                    "no G-related nuisance parameter estimates had to be bounded", 
                    # sum(apply(NuisanceDenom, 2, function(subj) any(subj < MinNuisance))), 
                    # "/", ncol(NuisanceDenom), " subjects", 
                    " to ", signif(MinNuisance, 3), sep = "")
            }
            
        }
    }
    return(NuisanceDenom)
}

