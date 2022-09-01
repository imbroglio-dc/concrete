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
    Censored <- 0 %in% Data[[attr(Data, "EventType")]]
    CovDT <- subset(Data, select = attr(Data, "CovNames")[["ColName"]])
    TrtModel <- try(Model[[attr(Data, "Treatment")]])
    if (inherits(TrtModel, "try-error"))
        stop("TrtModel must currently be specified in the Model argument as a list named ", 
             "as the Treatment variable, ", attr(Data, "Treatment"))
    
    ## Propensity Scores for Regimes of Interest ----
    PropScores <- getPropScore(TrtVal = TrtVal, CovDT = CovDT, TrtModel = TrtModel,
                               MinNuisance = MinNuisance, Regime = Regime,
                               PropScoreBackend = PropScoreBackend, CVFolds = CVFolds, TrtLoss = NULL, 
                               ReturnModels = ReturnModels)
    
    ## hazards: Events & censoring ----
    ## baseline hazards for obs times + target times ----
    HazTimes <- sort(unique(c(TargetTime, TimeVal)))
    HazTimes <- HazTimes[HazTimes <= max(TargetTime)]
    Hazards <- data.table("Time" = c(0, HazTimes))
    
    HazFits <- getHazFit(Data = Data, Model = Model, CVFolds = CVFolds, Hazards = Hazards, 
                         HazEstBackend = HazEstBackend, ReturnModels = ReturnModels)
    HazSurvPreds <- getHazSurvPred(Data = Data, HazFits = HazFits, MinNuisance = MinNuisance,
                                   TargetEvent = TargetEvent, TargetTime = TargetTime,
                                   Regime = Regime, HazEstBackend)
    InitialEstimates <- lapply(seq_along(PropScores), function(a) {
        if (Censored) {
            NuisanceWeight <- sapply(seq_along(PropScores[[a]]), function(i) {
                PropScores[[a]][i] * HazSurvPreds[[a]][["Survival"]][["LaggedCensSurv"]][, i]
            })   
        } else {
            NuisanceWeight <- sapply(seq_along(PropScores[[a]]), function(i) {
                PropScores[[a]][i] * 1
            })   
        }
        NuisanceWeight <- 1 / truncNuisanceDenom(NuisanceWeight, MinNuisance)
        return(list("PropScore" = PropScores[[a]],
                    "Hazards" = HazSurvPreds[[a]][["Hazards"]],
                    "EvntFreeSurv" = HazSurvPreds[[a]][["Survival"]][["TotalSurv"]],
                    "NuisanceWeight" = NuisanceWeight))
    })
    
    names(InitialEstimates) <- names(Regime)
    attr(InitialEstimates, "times") <- Hazards[["Time"]]
    if (ReturnModels) {
        ModelFits <- c(list(attr(PropScores, "TrtFit")), 
                       lapply(HazFits, function(HF) return(attr(HF, "HazSL"))))
        names(ModelFits)[1] <- attr(Data, "Treatment")
        attr(InitialEstimates, "ModelFits") <- ModelFits
    }
    return(InitialEstimates)
}

