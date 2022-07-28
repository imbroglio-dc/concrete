
#' Title
#'
#' @param Treatment numeric vector
#' @param CovDataTable data.table
#' @param TrtModel list or fitted object
#' @param MinNuisance numeric
#' @param Regime list
#' @param PropScoreBackend character
#' @param CVFolds list
#' @param TrtLoss character or function(A, g.A)
#' 
#' @import SuperLearner
#' @importFrom stats binomial gaussian

getPropScore <- function(Treatment, CovDataTable, TrtModel, MinNuisance, Regime,
                         PropScoreBackend, CVFolds, TrtLoss = NULL) {
    if (PropScoreBackend == "sl3") {
        PropScores <- getSl3PropScore(Treatment = Treatment, CovDataTable = CovDataTable,
                                      TrtModel = TrtModel, Regime = Regime, CVFolds = CVFolds)
    } else if (PropScoreBackend == "SuperLearner") {
        PropScores <- getSuperLearnerPropScore(Treatment = Treatment, CovDataTable = CovDataTable,
                                               TrtModel = TrtModel, Regime = Regime, cv.folds = CVFolds)
    } else {
        stop("functionality for propensity score calculation not using sl3 has not yet been implemented")
    }
    return(PropScores)
}

getSl3PropScore <- function(Treatment, CovDataTable, TrtModel, Regime, CVFolds) {
    ## PropScore score ----
    TrtTask <- sl3::make_sl3_Task(
        data = as.data.frame(cbind("Trt" = Treatment, CovDataTable)),
        covariates = colnames(CovDataTable),
        outcome = "Trt"
    )
    if (is.null(TrtModel$params$learners)) {
        TrtSL <- sl3::Lrnr_cv$new(learner = TrtModel, folds = CVFolds)
    } else {
        TrtSL <- sl3::Lrnr_sl$new(learners = TrtModel, folds = CVFolds)
    }
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
    attr(PropScores, "TrtFit") <- TrtFit
    return(PropScores)
}

getSuperLearnerPropScore <- function(Treatment, CovDataTable, TrtModel, Regime, cv.folds) {
    if (!inherits(TrtModel, "SuperLearner")) { 
        sl.args <- list() 
        sl.args[["Y"]] <- Treatment
        sl.args[["X"]] <- CovDataTable
        if (length(unique(Treatment) == 2))
            sl.args[["family"]] <- binomial()
        else
            sl.args[["family"]] <- gaussian()
        sl.args[["SL.library"]] <- TrtModel
        sl.args[["cvControl"]] <- list("V" = as.integer(length(cv.folds)), "stratifyCV" = FALSE, "shuffle" = FALSE,
                                       "validRows" = lapply(cv.folds, function(v) v[["validation_set"]]))
        superlearner.fit <- do.call(SuperLearner, sl.args)
    } else {
        superlearner.fit <- TrtModel
    }
    PropScores <- lapply(Regime, function(a) {
        if (is.numeric(a) & length(a) == length(Treatment)) {
            if (all(a %in% c(0, 1))) {
                ga1 <- as.numeric(predict.SuperLearner(object = superlearner.fit, newdata = CovDataTable)$pred)
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
