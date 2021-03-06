#' Title
#'
#' @param Data data.table
#' @param CovDataTable data.table
#' @param Model list
#' @param CVFolds : list
#' @param MinNuisance numeric
#' @param TargetEvent numeric vector
#' @param TargetTime numeric vector
#' @param Regime list
#' @param PropScoreBackend character
#' @param HazEstBackend character
#' @param Censored boolean

getInitialEstimate <- function(Data, CovDataTable, Model, CVFolds, MinNuisance, TargetEvent,
                               TargetTime, Regime, PropScoreBackend, HazEstBackend, Censored) {
    Time <- NULL
    Treatment <- Data[[attr(Data, "Treatment")]]
    EventType <- Data[[attr(Data, "EventType")]]
    EventTime <- Data[[attr(Data, "EventTime")]]

    TrtModel <- try(Model[[attr(Data, "Treatment")]])
    if (inherits(TrtModel, "try-error"))
        stop("TrtModel must currently be specified in the Model argument as a list named ", 
             "as the Treatment variable, ", attr(Data, "Treatment"))
    
    ## Propensity Scores for Regimes of Interest ----
    PropScores <- getPropScore(Treatment = Treatment, CovDataTable = CovDataTable, TrtModel = TrtModel,
                               MinNuisance = MinNuisance, Regime = Regime,
                               PropScoreBackend = PropScoreBackend, CVFolds = CVFolds, TrtLoss = NULL)

    ## hazards: Events & censoring ----
    ## baseline hazards for obs times + target times ----
    HazTimes <- sort(unique(c(TargetTime, EventTime)))
    HazTimes <- HazTimes[HazTimes <= max(TargetTime)]
    Hazards <- data.table("Time" = c(0, HazTimes))

    HazFits <- getHazFit(Data = Data, Model = Model, CVFolds = CVFolds, Hazards = Hazards, HazEstBackend)
    HazSurvPreds <- getHazSurvPred(Data = Data, HazFits = HazFits, MinNuisance = MinNuisance,
                                   TargetEvent = TargetEvent, TargetTime = TargetTime,
                                   Regime = Regime, HazEstBackend, Censored = Censored)
    InitialEstimates <- lapply(seq_along(PropScores), function(a) {
        NuisanceWeight <- sapply(seq_along(PropScores[[a]]), function(i) {
            PropScores[[a]][i] * HazSurvPreds[[a]][["Survival"]][["LaggedCensSurv"]][, i]
        })
        NuisanceWeight <- 1 / truncNuisanceDenom(NuisanceWeight, MinNuisance)
        return(list("PropScore" = PropScores[[a]],
                    "Hazards" = HazSurvPreds[[a]][["Hazards"]],
                    "EvntFreeSurv" = HazSurvPreds[[a]][["Survival"]][["TotalSurv"]],
                    "NuisanceWeight" = NuisanceWeight))
    })

    names(InitialEstimates) <- names(Regime)
    attr(InitialEstimates, "times") <- Hazards[["Time"]]
    return(InitialEstimates)
}

