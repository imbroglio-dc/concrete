#' Title
#'
#' @param Data data.table
#' @param CovDataTable data.table
#' @param Model list
#' @param MinNuisance numeric
#' @param TargetEvent numeric vector
#' @param TargetTime numeric vector
#' @param Regime list
#' @param PropScoreBackend character
#' @param Censored boolean
#'

getInitialEstimate <- function(Data, EventTime, EventType, Treatment, CovDataTable, Model, MinNuisance,
                               TargetEvent, TargetTime, Regime, PropScoreBackend, Censored) {
    Time <- NULL
    ## cross validation setup ----
    # stratifying cv so that folds are balanced for treatment assignment & outcomes
    # theory? but regressions may fail in practice with rare events otherwise ### make efficient CV representation ----
    StrataIDs <- as.numeric(factor(paste0(Treatment, ":", EventType)))
    CVFolds <- origami::make_folds(n = Data, fold_fun = origami::folds_vfold,
                                   strata_ids = StrataIDs)

    ## Propensity Scores for Regimes of Interest ----
    PropScores <- getPropScore(Treatment = Treatment,
                               CovDataTable, Model, MinNuisance, Regime,
                               PropScoreBackend, CVFolds)
    TrtFit <- PropScores[["TrtFit"]]
    PropScores <- PropScores[["PropScores"]]

    ## hazards: Events & censoring ----
    ## baseline hazards for obs times + target times ----
    HazTimes <- sort(unique(c(TargetTime, EventTime)))
    HazTimes <- HazTimes[HazTimes <= max(TargetTime)]
    Hazards <- data.table("Time" = c(0, HazTimes))

    HazFits <- getHazFit(Data = Data, EventTime = EventTime, Model = Model,
                         CVFolds = CVFolds, Hazards = Hazards)
    HazSurvPreds <- getHazSurvPred(Data, HazFits, MinNuisance, TargetEvent,
                                   TargetTime, Regime, Censored)
    InitialEstimates <- lapply(seq_along(PropScores), function(a) {
        NuisanceWeight <- sapply(seq_along(PropScores[[a]]), function(i) {
            PropScores[[a]][i] * HazSurvPreds[[a]][["Survival"]][["LaggedCensSurv"]][, i]})
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

