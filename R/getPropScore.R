
#' Title
#'
#' @param Treatment numeric vector
#' @param CovDataTable data.table
#' @param Model list
#' @param MinNuisance numeric
#' @param Regime list
#' @param PropScoreBackend character
#' @param CVFolds list
#' @param TrtLoss character or function(A, g.A)
#'

getPropScore <- function(Treatment, CovDataTable, Model, MinNuisance, Regime,
                         PropScoreBackend, CVFolds, TrtLoss = NULL) {
    if (PropScoreBackend == "sl3") {
        PropScores <- getSl3PropScore(Treatment, CovDataTable, Model, MinNuisance, Regime, CVFolds)
    } else {
        stop("functionality for propensity score calculation not using sl3 has not yet been implemented")
    }
    return(PropScores)
}

getSl3PropScore <- function(Treatment, CovDataTable, Model, MinNuisance, Regime, CVFolds) {
    ## PropScore score ----
    TrtTask <- sl3::make_sl3_Task(
        data = cbind("Trt" = Treatment, CovDataTable),
        covariates = colnames(CovDataTable),
        outcome = "Trt"
    )
    TrtSL <- sl3::Lrnr_sl$new(learners = Model[["Trt"]], folds = CVFolds)
    TrtFit <- TrtSL$train(TrtTask)

    PropScores <- lapply(Regime, function(a) {
        if (is.numeric(a) & length(a) == length(Treatment)) {
            if (all(a %in% c(0, 1))) {
                ga1 <- TrtFit$predict()
                PropScore <- ga1
                PropScore[a == 0] <- 1 - PropScore[a == 0]
            } else {
                stop("support for non-binary intervention variables is not yet implemented")
            }
        } else {
            stop("support for non-numeric, non-vector regimes of interest is not yet implemented.")
        }
        attr(PropScore, "g.star.intervention") <- attr(a, "g.star")(a)
        attr(PropScore, "g.star.obs") <- attr(a, "g.star")(Treatment)
        return(PropScore)
    })
    return(list("PropScores" = PropScores, "TrtFit" = TrtFit))
}

truncNuisanceDenom <- function(NuisanceDenom, MinNuisance) {
    if (is.function(MinNuisance))
        warning("Functionality for applying a MinNuisance function is not yet implemented")
    if (is.numeric(MinNuisance) & length(MinNuisance) == 1) {
        if (MinNuisance < 1 & MinNuisance > 0) {
            if (min(NuisanceDenom) < MinNuisance | max(NuisanceDenom) > 1 - MinNuisance) {
                warning("practical near positivity violation, truncating NuisanceDenom +-", MinNuisance, "\n")
                attr(NuisanceDenom, "original") <- NuisanceDenom
                NuisanceDenom[NuisanceDenom < MinNuisance] <- MinNuisance
                NuisanceDenom[NuisanceDenom > (1 - MinNuisance)] <- 1 - MinNuisance
            }
        }
    }
    return(NuisanceDenom)
}
