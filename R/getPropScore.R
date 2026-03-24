#' Estimate Propensity Scores for Treatment Assignment
#'
#' Estimates P(A|W), the conditional probability of treatment given covariates,
#' using cross-validated super learner. Computes both the propensity for the
#' observed treatment and for the intervention assignment (g.star).
#'
#' @param TrtVal data.table; observed treatment values (n x p for p treatment variables)
#' @param CovDT data.table; baseline covariates for prediction
#' @param TrtModel list; SuperLearner library specification for each treatment variable,
#'   or a pre-fitted SuperLearner object
#' @param MinNuisance numeric in (0, 1); minimum propensity score (for documentation;
#'   actual truncation occurs in getInitialEstimate)
#' @param Regime list; intervention specifications, each with g.star function
#' @param CVFolds list; cross-validation fold assignments
#' @param TrtLoss character or function (optional); loss function for super learner.
#'   Currently unused; defaults to binomial deviance for binary, MSE for continuous.
#' @param ReturnModels logical; if TRUE, store fitted SuperLearner objects
#'
#' @return List of propensity scores, one element per intervention regime.
#'   Each element is a numeric vector of length n with attributes:
#'   \itemize{
#'     \item g.star.intervention: P(A = a*|W) under the intervention
#'     \item g.star.obs: P(A = A_obs|W) for observed treatment
#'   }
#'   List attributes: TrtFit (fitted models if ReturnModels), warnings (any SuperLearner warnings)
#'
#' @details
#' For binary treatment: uses binomial family SuperLearner
#' For continuous treatment: uses Gaussian family (not fully implemented)
#' For multiple treatment variables: estimates joint propensity as product of conditionals
#'
#' The g.star function from each Regime computes the intervention-specific
#' propensity: typically 1(A = a*) for static interventions, but can be
#' more complex for dynamic or stochastic interventions.
#'
#' @import SuperLearner
#' @importFrom stats binomial gaussian
#' @keywords internal

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
                SLArgs[["family"]] <- ifelse(length(unique(unlist(TrtVal[[a_i]]))) == 2, "binomial", "gaussian")
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
                stop("Propensity score estimation requires either SuperLearner or sl3 backend. ",
                     "Received model of class: ", paste(class(TrtModel[[a_i]]), collapse = ", "), ". ",
                     "Use a character vector of SuperLearner library names (e.g., c('SL.glm', 'SL.glmnet')) ",
                     "or an sl3 learner object.",
                     call. = FALSE)
            }
        } 
    }
    
    PropScores <- lapply(Regime, function(a) {
        if (!all(dim(a) == dim(TrtVal))) {
            stop("Intervention regime dimensions (", paste(dim(a), collapse = " x "),
                 ") do not match observed treatment dimensions (",
                 paste(dim(TrtVal), collapse = " x "), "). ",
                 "The intervention function must return a data structure with the same ",
                 "dimensions as the treatment variable. Check your Intervention specification.",
                 call. = FALSE)
        }
        
        PropScore <- rep_len(1, nrow(TrtVal))
        for (a_i in 1:ncol(TrtVal)) {
            a_vec <- unlist(subset(a, select = a_i))
            if (!all(a_vec %in% c(0, 1))) {
                stop("Non-binary intervention variables are not yet supported. ",
                     "Intervention values must be 0 or 1. Found values: ",
                     paste(unique(a_vec[!a_vec %in% c(0, 1)])[1:min(5, sum(!a_vec %in% c(0, 1)))], collapse = ", "),
                     ". For continuous or multi-level treatments, consider dichotomizing ",
                     "or using a different estimation approach.",
                     call. = FALSE)
            }
            
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

