#' Title
#'
#' @param Data data.table
#' @param CovDataTable data.table
#' @param Models list
#' @param MinNuisance numeric
#' @param TargetEvents numeric vector
#' @param TargetTimes numeric vector
#' @param RegsOfInterest list
#' @param PropScoreBackend character
#' @param Censored boolean
#'

getInitialEstimate <- function(Data, CovDataTable, Models, MinNuisance, TargetEvents,
                               TargetTimes, RegsOfInterest, PropScoreBackend, Censored) {
    Time <- NULL
    ## cross validation setup ----
    # stratifying cv so that folds are balanced for treatment assignment & outcomes
    # theory? but regressions may fail in practice with rare events otherwise
    StrataIDs <- as.numeric(factor(paste0(Data[["Trt"]], ":", Data[["Event"]])))
    CVFolds <- origami::make_folds(n = Data, fold_fun = origami::folds_vfold,
                                   strata_ids = StrataIDs)

    ## Propensity Scores for Regimes of Interest ----
    PropScores <- getPropScore(Treatment = Data[["Trt"]],
                               CovDataTable, Models, MinNuisance, RegsOfInterest,
                               PropScoreBackend, CVFolds)
    TrtFit <- PropScores[["TrtFit"]]
    PropScores <- PropScores[["PropScores"]]

    ## hazards: Events & censoring ----
    ## baseline hazards for obs times + target times ----
    HazTimes <- unique(c(TargetTimes, Data[["Time"]]))
    HazTimes <- HazTimes[HazTimes <= max(TargetTimes)]
    Hazards <- data.table("Time" = c(0, HazTimes))[order(Time)]

    HazFits <- getHazFit(Data, Models, CVFolds, Hazards)
    HazSurvPreds <- getHazSurvPred(Data, HazFits, MinNuisance, TargetEvents,
                                   TargetTimes, RegsOfInterest, Censored)
    InitialEstimates <- lapply(seq_along(PropScores), function(a) {
        NuisanceWeight <- sapply(seq_along(PropScores[[a]]), function(i) {
            PropScores[[a]][i] * HazSurvPreds[[a]][["Survival"]][["LaggedCensSurv"]][, i]})
        NuisanceWeight <- 1 / truncNuisanceDenom(NuisanceWeight, MinNuisance)
        return(list("PropScore" = PropScores[[a]],
                    "Hazards" = HazSurvPreds[[a]][["Hazards"]],
                    "EvntFreeSurv" = HazSurvPreds[[a]][["Survival"]][["TotalSurv"]],
                    "NuisanceWeight" = NuisanceWeight))
    })

    names(InitialEstimates) <- names(RegsOfInterest)
    attr(InitialEstimates, "times") <- Hazards[["Time"]]
    return(InitialEstimates)
}

