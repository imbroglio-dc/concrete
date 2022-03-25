
#' Title
#'
#' @param Data data.table
#' @param CovDataTable data.table
#' @param Models list
#' @param MinNuisance numeric
#' @param RegsOfInterest list
#' @param PropScoreBackend character
#' @param CVFolds list
#'

getPropScores <- function(Data, CovDataTable, Models, MinNuisance, RegsOfInterest,
                          PropScoreBackend, CVFolds) {
    if (PropScoreBackend == "sl3") {
        PropScores <- getSl3PropScores(Data, CovDataTable, Models, MinNuisance, RegsOfInterest, CVFolds)
    } else {
        stop("functionality for propensity score calculation not using sl3 has not yet been implemented")
    }
    return(PropScores)
}

getSl3PropScores <- function(Data, CovDataTable, Models, MinNuisance, RegsOfInterest, CVFolds) {
    ## PropScore score ----
    TrtTask <- sl3::make_sl3_Task(
        data = Data[, -c("Time", "Event", "ID")],
        covariates = colnames(CovDataTable),
        outcome = "Trt"
    )
    TrtSL <- sl3::Lrnr_sl$new(learners = Models[["Trt"]], folds = CVFolds)
    TrtFit <- TrtSL$train(TrtTask)

    PropScores <- lapply(RegsOfInterest, function(a) {
        if (is.numeric(a) & length(a) == length(Data[["Trt"]])) {
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
        attr(PropScore, "g.star.obs") <- attr(a, "g.star")(Data[["Trt"]])
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
