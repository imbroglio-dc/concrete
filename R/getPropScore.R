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
#' @param Regime list; intervention specifications, each with a g.star attribute
#' @param CVFolds list; cross-validation fold assignments
#' @param TrtLoss character or function (optional); loss function for super learner.
#'   Currently unused; defaults to binomial deviance for binary, MSE for continuous.
#' @param ReturnModels logical; if TRUE, store fitted SuperLearner objects
#'
#' @return List of propensity scores, one element per intervention regime.
#'   Each element is a numeric vector of length n with attributes:
#'   \itemize{
#'     \item \code{g.star.obs}: g*(A_obs | W) - P(observed A matches regime | W).
#'       Used in the efficient influence curve (an observed-data score).
#'     \item \code{g.star.intervention}: g*(A* | W) - P(regime assignment | W) under
#'       the target intervention. Used for the TMLE hazard update direction, which
#'       fluctuates the counterfactual hazards of every subject. For static ITT
#'       interventions this is identically 1; for stochastic interventions it varies
#'       across subjects.
#'   }
#'   List attributes:
#'   \itemize{
#'     \item \code{TrtFit}: fitted SuperLearner object(s) (if \code{ReturnModels = TRUE}),
#'       or a descriptive string otherwise
#'     \item \code{SLWarnings}: named list (one element per treatment variable) of any
#'       warnings emitted during SuperLearner fitting; empty character vectors when none
#'   }
#'
#' @details
#' For binary treatment: uses binomial family SuperLearner.
#' For continuous treatment: uses Gaussian family (not fully implemented; use with caution).
#' For multiple treatment variables: estimates joint propensity as product of conditionals.
#'
#' Raw P(A=1|W) predictions are computed once per treatment variable and reused across
#' all intervention regimes, avoiding redundant SuperLearner predict calls.
#'
#' @import SuperLearner
#' @importFrom stats binomial gaussian
#' @keywords internal

getPropScore <- function(TrtVal, CovDT, TrtModel, MinNuisance, Regime,
                         CVFolds, TrtLoss = NULL, ReturnModels) {

    # --- guard: degenerate treatment variables --------------------------------
    for (a_i in seq_len(ncol(TrtVal))) {
        trt_name <- colnames(TrtVal)[a_i]
        n_unique  <- length(unique(unlist(TrtVal[[a_i]])))
        if (n_unique < 2L) {
            stop("Treatment column '", trt_name, "' has only one unique value; ",
                 "propensity score estimation requires at least 2 distinct treatment values. ",
                 "Check that treatment assignment varies across subjects.",
                 call. = FALSE)
        }
    }

    # --- fit (or reuse) propensity models ------------------------------------
    # SLWarnings collects warnings emitted during SuperLearner fitting, scoped
    # per treatment variable (not the entire session, unlike warnings()).
    SLWarnings <- vector("list", ncol(TrtVal))
    names(SLWarnings) <- colnames(TrtVal)

    if (all(sapply(TrtModel, function(a) inherits(a, "SuperLearner")))) {
        TrtFit <- TrtModel
    } else {
        TrtFit <- vector("list", ncol(TrtVal))
        names(TrtFit) <- colnames(TrtVal)

        for (a_i in seq_len(ncol(TrtVal))) {
            trt_name <- colnames(TrtVal)[a_i]
            if (attr(TrtModel[[a_i]], "Backend") == "SuperLearner") {
                SLArgs <- list()
                SLArgs[["Y"]] <- unlist(subset(TrtVal, select = a_i))
                SLArgs[["X"]] <- if (a_i > 1L) {
                    cbind(subset(TrtVal, select = seq_len(a_i - 1L)), CovDT)
                } else {
                    CovDT
                }
                n_unique_y <- length(unique(SLArgs[["Y"]]))
                if (n_unique_y == 2L) {
                    SLArgs[["family"]] <- "binomial"
                } else {
                    SLArgs[["family"]] <- "gaussian"
                    warning("Treatment '", trt_name, "' has ", n_unique_y,
                            " unique values; fitting propensity score with Gaussian family. ",
                            "Continuous treatment support is incomplete - interpret results carefully.",
                            call. = FALSE)
                }
                # Strip non-standard attributes (e.g., "Backend") before passing to SuperLearner
                SLArgs[["SL.library"]] <- as.character(TrtModel[[a_i]])
                SLArgs[["cvControl"]] <- list("V"          = as.integer(length(CVFolds)),
                                              "stratifyCV" = FALSE,
                                              "shuffle"    = FALSE,
                                              "validRows"  = lapply(CVFolds, `[[`, "validation_set"))

                fit_warnings <- character(0L)
                TrtFit[[a_i]] <- tryCatch(
                    withCallingHandlers(
                        do.call(SuperLearner, SLArgs),
                        warning = function(w) {
                            fit_warnings <<- c(fit_warnings, conditionMessage(w))
                            invokeRestart("muffleWarning")
                        }
                    ),
                    error = function(e) {
                        stop("SuperLearner propensity score estimation failed for treatment '",
                             trt_name, "' (variable ", a_i, " of ", ncol(TrtVal), "). ",
                             "SL library: ", paste(SLArgs[["SL.library"]], collapse = ", "), ". ",
                             "Original error: ", conditionMessage(e), "\n",
                             "Consider simplifying the SL library or checking the CV fold structure.",
                             call. = FALSE)
                    }
                )
                SLWarnings[[trt_name]] <- fit_warnings
            # } else if (attr(TrtModel[[a_i]], "Backend") == "sl3") { ... }
            } else {
                stop("Propensity score estimation requires either SuperLearner or sl3 backend. ",
                     "Received model of class: ", paste(class(TrtModel[[a_i]]), collapse = ", "), ". ",
                     "Use a character vector of SuperLearner library names (e.g., c('SL.glm', 'SL.glmnet')) ",
                     "or an sl3 learner object.",
                     call. = FALSE)
            }
        }
    }

    # --- precompute raw P(A=1|W) once, outside the regime loop ---------------
    # predict.SuperLearner is called once per treatment variable regardless of
    # how many regimes there are; the per-regime flip is applied below.
    g_raw <- vector("list", ncol(TrtVal))
    names(g_raw) <- colnames(TrtVal)

    for (a_i in seq_len(ncol(TrtVal))) {
        trt_name <- colnames(TrtVal)[a_i]
        if (!inherits(TrtFit[[a_i]], "SuperLearner")) {
            stop("TrtFit[[", a_i, "]] ('", trt_name, "') is not a SuperLearner object ",
                 "(class: ", paste(class(TrtFit[[a_i]]), collapse = ", "), "). ",
                 "This is unexpected; please file a bug report.",
                 call. = FALSE)
        }
        newdata <- if (a_i > 1L) {
            cbind(subset(TrtVal, select = seq_len(a_i - 1L)), CovDT)
        } else {
            CovDT
        }
        g_pred <- tryCatch(
            as.numeric(predict.SuperLearner(object = TrtFit[[a_i]], newdata = newdata)$pred),
            error = function(e) {
                stop("SuperLearner prediction failed for treatment '", trt_name, "'. ",
                     "Original error: ", conditionMessage(e),
                     call. = FALSE)
            }
        )
        if (anyNA(g_pred) || any(is.nan(g_pred))) {
            n_bad <- sum(is.na(g_pred) | is.nan(g_pred))
            stop("SuperLearner produced ", n_bad, " NA/NaN prediction(s) for treatment '",
                 trt_name, "'. This usually indicates a learner failure within the ensemble. ",
                 "Try simplifying the SL library or increasing MinNuisance.",
                 call. = FALSE)
        }
        g_raw[[a_i]] <- g_pred
    }

    # --- compute per-regime propensity scores --------------------------------
    PropScores <- lapply(seq_along(Regime), function(a_regime) {
        a           <- Regime[[a_regime]]
        regime_name <- names(Regime)[a_regime]

        if (!all(dim(a) == dim(TrtVal))) {
            stop("Intervention regime '", regime_name, "' dimensions (",
                 paste(dim(a), collapse = " x "),
                 ") do not match observed treatment dimensions (",
                 paste(dim(TrtVal), collapse = " x "), "). ",
                 "The intervention function must return a data structure with the same ",
                 "dimensions as the treatment variable. Check your Intervention specification.",
                 call. = FALSE)
        }

        PropScore <- rep_len(1, nrow(TrtVal))
        for (a_i in seq_len(ncol(TrtVal))) {
            a_vec <- unlist(subset(a, select = a_i))
            if (!all(a_vec %in% c(0, 1))) {
                bad_vals <- unique(a_vec[!a_vec %in% c(0, 1)])
                stop("Non-binary intervention values in regime '", regime_name, "'. ",
                     "Intervention values must be 0 or 1. Found: ",
                     paste(bad_vals[seq_len(min(5L, length(bad_vals)))], collapse = ", "),
                     ". For continuous or multi-level treatments, consider dichotomizing ",
                     "or using a different estimation approach.",
                     call. = FALSE)
            }
            # Use precomputed raw P(A=1|W); flip probability for subjects assigned to control
            g.a          <- g_raw[[a_i]]
            g.a[a_vec == 0L] <- 1 - g.a[a_vec == 0L]
            PropScore    <- PropScore * g.a
        }

        g_star_fn <- attr(a, "g.star")
        if (!is.function(g_star_fn)) {
            stop("Regime '", regime_name, "' has no g.star function. ",
                 "Ensure the regime was created via makeITT() or includes a valid g.star attribute.",
                 call. = FALSE)
        }

        # g.star.obs: g*(A_obs | W) - used in the EIC (observed-data score)
        attr(PropScore, "g.star.obs") <- tryCatch(
            g_star_fn(TrtVal, CovDT, PropScore, a),
            error = function(e) {
                stop("g.star function failed for regime '", regime_name, "' ",
                     "(observed treatment). Original error: ", conditionMessage(e),
                     call. = FALSE)
            }
        )

        # g.star.intervention: g*(A* | W) - probability that the intervention assigns
        # the target treatment value to each subject. Used for the TMLE hazard update
        # direction. For static ITT interventions this is identically 1; for
        # stochastic interventions it varies across subjects.
        attr(PropScore, "g.star.intervention") <- tryCatch(
            g_star_fn(a, CovDT, PropScore, a),
            error = function(e) {
                stop("g.star function failed for regime '", regime_name, "' ",
                     "(intervention treatment). Original error: ", conditionMessage(e),
                     call. = FALSE)
            }
        )

        return(PropScore)
    })
    names(PropScores) <- names(Regime)

    attr(PropScores, "TrtFit")    <- if (ReturnModels) TrtFit else
        "TrtFits not saved because `ReturnModels' was set = FALSE"
    attr(PropScores, "SLWarnings") <- Filter(length, SLWarnings)
    return(PropScores)
}
