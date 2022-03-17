
getPropScores <- function(Data, CovDataTable, Models, MinNuisanceDenom, RegsOfInterest,
                          PropScoreBackend, CVFolds) {
    if (PropScoreBackend == "sl3") {
        PropScores <- getSl3PropScores(Data, CovDataTable, Models, MinNuisanceDenom, RegsOfInterest, CVFolds)
    } else {
        stop("functionality for propensity score calculation not using sl3 has not yet been implemented")
    }
}

getSl3PropScores <- function(Data, CovDataTable, Models, MinNuisanceDenom, RegsOfInterest, CVFolds) {
    ## PropScore score ----
    TrtTask <- sl3::make_sl3_Task(
        data = Data[, -c("Time", "Event", "ID")],
        covariates = colnames(CovDataTable),
        outcome = "Trt"
    )

    TrtSL <- sl3::Lrnr_sl$new(learners = Models[["A"]], folds = CVFolds)
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
        return(PropScore)
    })
    return(PropScores)
}

truncNuisanceDenom <- function(NuisanceDenom, MinNuisanceDenom) {
    if (is.function(MinNuisanceDenom))
        warning("Functionality for applying a MinNuisanceDenom function is not yet implemented")
    if (is.numeric(MinNuisanceDenom) & length(MinNuisanceDenom) == 1) {
        if (MinNuisanceDenom < 1 & MinNuisanceDenom > 0) {
            if (min(NuisanceDenom) < MinNuisanceDenom | max(NuisanceDenom) > 1 - MinNuisanceDenom) {
                warning("practical near positivity violation, truncating NuisanceDenom +-", MinNuisanceDenom, "\n")
                attr(NuisanceDenom, "original") <- NuisanceDenom
                NuisanceDenom[NuisanceDenom < MinNuisanceDenom] <- MinNuisanceDenom
                NuisanceDenom[NuisanceDenom > (1 - MinNuisanceDenom)] <- 1 - MinNuisanceDenom
            }
        }
    }
    return(NuisanceDenom)
}
