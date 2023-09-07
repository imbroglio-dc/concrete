#' getPropScore
#'
#' @param TrtVal numeric vector
#' @param CovDT data.table
#' @param TrtModel list or fitted object
#' @param MinNuisance numeric
#' @param Regime list
#' @param CVFolds list
#' @param TrtLoss character or function(A, g.A)
#' @param ReturnModels boolean
#' 
#' @import SuperLearner
#' @importFrom stats binomial gaussian

getPropScore <- function(TrtVal, CovDT, TrtModel, MinNuisance, Regime, 
                         CVFolds, TrtLoss = NULL, ReturnModels) {
    old <- options()
    on.exit(options(old))
    options(warn = 0)
    
    if (all(sapply(TrtModel, function(a) inherits(a, "SuperLearner")))) {
        TrtFit <- TrtModel
    } else {
        TrtFit <- vector("list", ncol(TrtVal))
        names(TrtFit) <- colnames(TrtVal)
        for (a_i in 1:ncol(TrtVal)) {
            if (attr(TrtModel[[a_i]], "Backend") == "SuperLearner") {
                SLArgs <- list()
                SLArgs[["Y"]] <- unlist(subset(TrtVal, select = a_i))
                if (a_i > 1) {
                    SLArgs[["X"]] <- cbind(subset(TrtVal, select = 1:(a_i - 1)), CovDT)
                } else {
                    SLArgs[["X"]] <- CovDT
                }
                SLArgs[["family"]] <- ifelse(length(unique(unlist(TrtVal))) == 2, "binomial", "gaussian")
                SLArgs[["SL.library"]] <- TrtModel[[a_i]]
                SLArgs[["cvControl"]] <- list("V" = as.integer(length(CVFolds)), "stratifyCV" = FALSE, 
                                              "shuffle" = FALSE,
                                              "validRows" = lapply(CVFolds, function(v) v[["validation_set"]]))
                TrtFit[[a_i]] <- do.call(SuperLearner, SLArgs)
            # } else if (attr(TrtModel[[a_i]], "Backend") == "sl3") {
            #     data <- as.data.frame(cbind(subset(TrtVal, select = 1:a_i), CovDT))
            #     TrtTask <- sl3::make_sl3_Task(
            #         data = data,
            #         covariates = setdiff(colnames(data), colnames(TrtVal)[a_i]),
            #         outcome = colnames(TrtVal)[a_i]
            #     )
            #     if (is.null(TrtModel[[a_i]]$params$learners)) {
            #         TrtSL <- sl3::Lrnr_cv$new(learner = TrtModel[[a_i]], folds = CVFolds)
            #     } else {
            #         TrtSL <- sl3::Lrnr_sl$new(learners = TrtModel[[a_i]], folds = CVFolds)
            #     }
            #     TrtFit[[a_i]] <- TrtSL$train(TrtTask)
            } else {
                stop("functionality for propensity score estimation not using 'sl3' or ", 
                     "'SuperLearner' has not yet been implemented")
            }
        } 
    }
    
    PropScores <- lapply(Regime, function(a) {
        if (!all(dim(a) == dim(TrtVal))) 
            stop("Regime dimensions don't match with observed treatment. Bugfix needed")
        
        PropScore <- rep_len(1, nrow(TrtVal))
        for (a_i in 1:ncol(TrtVal)) {
            a_vec <- unlist(subset(a, select = a_i))
            if (!all(a_vec %in% c(0, 1))) 
                stop("support for non-binary intervention variables is not yet implemented")
            
            if (attr(TrtModel[[a_i]], "Backend") == "SuperLearner") {
                if (a_i > 1) {
                    newdata <- cbind(subset(TrtVal, select = 1:(a_i - 1)), CovDT)
                } else {
                    newdata <- CovDT
                }
                g.a <- as.numeric(predict.SuperLearner(object = TrtFit[[a_i]], newdata = newdata)$pred)
                g.a[a_vec == 0] <- 1 - g.a[a_vec == 0]
                PropScore <- PropScore * g.a
            # } else if (attr(TrtModel[[a_i]], "Backend") == "sl3") {
            #     g.a <- unlist(TrtFit[[a_i]]$predict())
            #     g.a[a_vec == 0] <- 1 - g.a[a_vec == 0]
            #     PropScore <- PropScore * g.a
            } 
        }
        attr(PropScore, "g.star.intervention") <- attr(a, "g.star")(a, CovDT, PropScore, a)
        attr(PropScore, "g.star.obs") <- attr(a, "g.star")(TrtVal, CovDT, PropScore, a)
        return(PropScore)
    })
    
    if (ReturnModels) {
        attr(PropScores, "TrtFit") <- TrtFit
    } else {
        attr(PropScores, "TrtFit") <- "TrtFits not saved because `ReturnModels' was set = FALSE"
    }
    attr(PropScores, "warnings") <- summary(warnings())
    return(PropScores)
}

