
#' Title
#'
#' @param TrtVal numeric vector
#' @param CovDT data.table
#' @param TrtModel list or fitted object
#' @param MinNuisance numeric
#' @param Regime list
#' @param PropScoreBackend character
#' @param CVFolds list
#' @param TrtLoss character or function(A, g.A)
#' @param ReturnModels boolean
#' 
#' @import SuperLearner
#' @importFrom stats binomial gaussian

getPropScore <- function(TrtVal, CovDT, TrtModel, MinNuisance, Regime,
                         PropScoreBackend, CVFolds, TrtLoss = NULL, ReturnModels) {
    if (PropScoreBackend == "sl3") {
        PropScores <- getSl3PropScore(TrtVal = TrtVal, CovDT = CovDT, TrtModel = TrtModel, 
                                      Regime = Regime, CVFolds = CVFolds, ReturnModels = ReturnModels)
    } else if (PropScoreBackend == "SuperLearner") {
        PropScores <- getSuperLearnerPropScore(TrtVal = TrtVal, CovDT = CovDT, TrtModel = TrtModel, 
                                               Regime = Regime, CVFolds = CVFolds, ReturnModels = ReturnModels)
    } else {
        stop("functionality for propensity score calculation not using sl3 has not yet been implemented")
    }
    return(PropScores)
}

getSl3PropScore <- function(TrtVal, CovDT, TrtModel, Regime, CVFolds, ReturnModels) {
    ## PropScore score ----
    TrtTask <- sl3::make_sl3_Task(
        data = as.data.frame(cbind("Trt" = TrtVal, CovDT)),
        covariates = colnames(CovDT),
        outcome = "Trt"
    )
    if (is.null(TrtModel$params$learners)) {
        TrtSL <- sl3::Lrnr_cv$new(learner = TrtModel, folds = CVFolds)
    } else {
        TrtSL <- sl3::Lrnr_sl$new(learners = TrtModel, folds = CVFolds)
    }
    TrtFit <- TrtSL$train(TrtTask)
    
    PropScores <- lapply(Regime, function(a) {
        if (is.numeric(a) & length(a) == length(TrtVal)) {
            if (all(a %in% c(0, 1))) {
                ga1 <- unlist(TrtFit$predict())
                PropScore <- ga1
                PropScore[a == 0] <- 1 - PropScore[a == 0]
            } else {
                stop("support for non-binary intervention variables is not yet implemented")
            }
        } else {
            stop("support for non-numeric, non-vector regimes of interest is not yet implemented.")
        }
        attr(PropScore, "g.star.intervention") <- attr(a, "g.star")(a, CovDT)
        attr(PropScore, "g.star.obs") <- attr(a, "g.star")(TrtVal, CovDT)
        return(PropScore)
    })
    if (ReturnModels) attr(PropScores, "TrtFit") <- TrtFit
    return(PropScores)
}

getSuperLearnerPropScore <- function(TrtVal, CovDT, TrtModel, Regime, CVFolds, ReturnModels) {
    if (!inherits(TrtModel, "SuperLearner")) { 
        SLArgs <- list() 
        SLArgs[["Y"]] <- TrtVal
        SLArgs[["X"]] <- CovDT
        SLArgs[["family"]] <- gaussian()
        if (length(unique(TrtVal) == 2)) 
            SLArgs[["family"]] <- binomial()
        SLArgs[["SL.library"]] <- TrtModel
        SLArgs[["cvControl"]] <- list("V" = as.integer(length(CVFolds)), "stratifyCV" = FALSE, "shuffle" = FALSE,
                                       "validRows" = lapply(CVFolds, function(v) v[["validation_set"]]))
        TrtFit <- do.call(SuperLearner, SLArgs)
    } else {
        TrtFit <- TrtModel
    }
    PropScores <- lapply(Regime, function(a) {
        if (is.numeric(a) & length(a) == length(TrtVal)) {
            if (all(a %in% c(0, 1))) {
                ga1 <- as.numeric(predict.SuperLearner(object = TrtFit, newdata = CovDT)$pred)
                PropScore <- ga1
                PropScore[a == 0] <- 1 - PropScore[a == 0]
            } else {
                stop("support for non-binary intervention variables is not yet implemented")
            }
        } else {
            stop("support for non-numeric, non-vector regimes of interest is not yet implemented.")
        }
        attr(PropScore, "g.star.intervention") <- attr(a, "g.star")(a, CovDT)
        attr(PropScore, "g.star.obs") <- attr(a, "g.star")(TrtVal, CovDT)
        
        return(PropScore)
    })
    if (ReturnModels) attr(PropScores, "TrtFit") <- TrtFit
    return(PropScores)
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
