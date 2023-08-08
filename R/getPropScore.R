
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
    options(warn = 0)
    if (PropScoreBackend == "sl3") {
        PropScores <- getSl3PropScore(TrtVal = TrtVal, CovDT = CovDT, TrtModel = TrtModel, 
                                      Regime = Regime, CVFolds = CVFolds, ReturnModels = ReturnModels)
    } else if (PropScoreBackend == "SuperLearner") {
        PropScores <- getSuperLearnerPropScore(TrtVal = TrtVal, CovDT = CovDT, TrtModel = TrtModel, 
                                               Regime = Regime, CVFolds = CVFolds, ReturnModels = ReturnModels)
    } else {
        stop("functionality for propensity score calculation not using sl3 has not yet been implemented")
    }
    attr(PropScores, "warnings") <- summary(warnings())
    last.warning <- NULL
    return(PropScores)
}

getSl3PropScore <- function(TrtVal, CovDT, TrtModel, Regime, CVFolds, ReturnModels) {
    for (a_i in 1:ncol(TrtVal)) {
        ## PropScore score ----
        data <- as.data.frame(cbind(subset(TrtVal, select = 1:a_i), CovDT))
        TrtTask <- sl3::make_sl3_Task(
            data = data,
            covariates = setdiff(colnames(data), colnames(TrtVal)[a_i]),
            outcome = colnames(TrtVal)[a_i]
        )
        if (is.null(TrtModel$params$learners)) {
            TrtSL <- sl3::Lrnr_cv$new(learner = TrtModel, folds = CVFolds)
        } else {
            TrtSL <- sl3::Lrnr_sl$new(learners = TrtModel, folds = CVFolds)
        }
        TrtFit <- TrtSL$train(TrtTask)
    }
    
    
    PropScores <- lapply(Regime, function(a) {
        browser()
        if (is.numeric(a) & length(a) == nrow(TrtVal)) {
            for (a_i in 1:ncol(TrtVal)) {
                if (all(a %in% c(0, 1))) {
                    ga1 <- unlist(TrtFit$predict())
                    PropScore <- ga1
                    PropScore[a == 0] <- 1 - PropScore[a == 0]
                } else {
                    stop("support for non-binary intervention variables is not yet implemented")
                }
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
    browser()
    if (inherits(TrtModel, "SuperLearner")) {
        TrtFit <- TrtModel
    } else {
        TrtFit <- list()
        for (a_i in 1:ncol(TrtVal)) {
            SLArgs <- list()
            SLArgs[["Y"]] <- unlist(subset(TrtVal, select = a_i))
            if (a_i > 1) {
                SLArgs[["X"]] <- cbind(subset(TrtVal, select = 1:(a_i - 1)), CovDT)
            } else {
                SLArgs[["X"]] <- CovDT
            }
            SLArgs[["family"]] <- ifelse(length(unique(TrtVal)) == 2, binomial(), gaussian())
            SLArgs[["SL.library"]] <- TrtModel
            SLArgs[["cvControl"]] <- list("V" = as.integer(length(CVFolds)), "stratifyCV" = FALSE, "shuffle" = FALSE,
                                          "validRows" = lapply(CVFolds, function(v) v[["validation_set"]]))
            TrtFit[[a_i]] <- do.call(SuperLearner, SLArgs)
        }
    }
    
    PropScores <- lapply(Regime, function(a) {
        if (dim(a) == dim(TrtVal)) { 
            for (a_i in 1:ncol(TrtVal)) {
                a_vec <- unlist(subset(Regime, select = a_i))
                if (all(a_vec %in% c(0, 1))) {
                    ga1 <- as.numeric(predict.SuperLearner(object = TrtFit, newdata = CovDT)$pred)
                    PropScore <- ga1
                    PropScore[a_vec == 0] <- 1 - PropScore[a_vec == 0]
                } else
                    stop("support for non-binary intervention variables is not yet implemented")
            }
        } else
            stop("Regime dimensions don't match with observed treatment. Something went badly wrong")
        
        attr(PropScore, "g.star.intervention") <- attr(a, "g.star")(a, CovDT)
        attr(PropScore, "g.star.obs") <- attr(a, "g.star")(TrtVal, CovDT)
        
        return(PropScore)
    })
    if (ReturnModels) 
        attr(PropScores, "TrtFit") <- TrtFit
    else 
        attr(PropScores, "TrtFit") <- cbind(Risk = TrtFit$cvRisk, Coef = TrtFit$coef)
    return(PropScores)
}
