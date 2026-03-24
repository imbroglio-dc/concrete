#' Compute Initial Estimates for TMLE
#'
#' Estimates the initial (pre-targeting) nuisance parameters required for TMLE:
#' propensity scores, cause-specific hazards, and overall survival. These estimates
#' are computed using cross-validated super learner ensembles.
#'
#' @param Data data.table; formatted data with column name attributes
#' @param Model list; model specifications for propensity and hazard estimation.
#'   Named elements for Treatment column(s) specify SuperLearner libraries;
#'   elements named by EventType values specify Cox model formulas or survivalSL specs.
#' @param CVFolds list; cross-validation fold assignments from origami::make_folds()
#' @param MinNuisance numeric in (0, 1); minimum value for propensity score truncation.
#'   Values below this threshold are set to MinNuisance to avoid numerical instability.
#' @param TargetEvent numeric vector; event types to target for risk estimation
#' @param TargetTime numeric vector; time points for risk estimation
#' @param Regime list; intervention specifications, each with a g.star function
#' @param ReturnModels logical; if TRUE, fitted model objects are stored for inspection
#'
#' @return List with one element per intervention regime, each containing:
#'   \itemize{
#'     \item PropScore: propensity score estimates P(A=a|W)
#'     \item Hazards: list of cause-specific hazard matrices (times x subjects)
#'     \item EvntFreeSurv: overall survival S(t) matrix
#'     \item NuisanceWeight: inverse of truncated g-related weights (for EIC)
#'   }
#'   With attributes: Times (evaluation time grid), InitFits (fitted models if ReturnModels=TRUE)
#'
#' @details
#' The nuisance weight combines propensity scores and censoring survival:
#' weight = g(A|W) * S_C(t|A,W), truncated below by MinNuisance to ensure
#' finite variance of the EIC. This truncation is the primary source of
#' potential bias-variance tradeoff in TMLE.
#'
#' @keywords internal
getInitialEstimate <- function(Data, Model, CVFolds, MinNuisance, TargetEvent, TargetTime,
                               Regime, ReturnModels) {
    Time <- NULL
    TimeVal <- Data[[attr(Data, "EventTime")]]
    Censored <- any(Data[[attr(Data, "EventType")]] <= 0)
    
    ## Propensity Scores for Regimes of Interest ----
    message("\nEstimating Treatment Propensity:\n")
    PropScores <- getPropScore(TrtVal = Data[, .SD, .SDcols = attr(Data, "Treatment")], 
                               CovDT = subset(Data, select = attr(Data, "CovNames")[["ColName"]]), 
                               TrtModel = Model[which(names(Model) %in% attr(Data, "Treatment"))],
                               MinNuisance = MinNuisance, Regime = Regime, CVFolds = CVFolds, 
                               TrtLoss = NULL, ReturnModels = ReturnModels)
    InitFits <- attr(PropScores, "TrtFit")
    
    ## hazards: Events & censoring ----
    ## baseline hazards for obs times + target times ----
    HazTimes <- sort(unique(c(TargetTime, TimeVal)))
    HazTimes <- HazTimes[HazTimes <= max(TargetTime)]
    Hazards <- data.table("Time" = c(0, HazTimes))
    
    message("\nEstimating Hazards:\n")
    HazFits <- getHazFit(Data = Data, Model = Model, CVFolds = CVFolds, Hazards = Hazards, 
                         ReturnModels = ReturnModels)
    InitFits <- c(InitFits, lapply(HazFits, function(HF) return(attr(HF, "HazSL"))))
    HazSurvPreds <- getHazSurvPred(Data = Data, HazFits = HazFits, MinNuisance = MinNuisance,
                                   TargetEvent = TargetEvent, TargetTime = TargetTime,
                                   Regime = Regime)
    InitialEstimates <- lapply(seq_along(PropScores), function(a) {
        if (Censored) {
            NuisanceDenom <- sapply(seq_along(PropScores[[a]]), function(i) {
                PropScores[[a]][i] * HazSurvPreds[[a]][["Survival"]][["LaggedCensSurv"]][, i]
            })   
        } else {
            Srv <- HazSurvPreds[[a]][["Survival"]][["TotalSurv"]]
            NuisanceDenom <- matrix(PropScores[[a]], nrow = nrow(Srv), ncol = ncol(Srv), byrow = TRUE) 
        }
        NuisanceWeight <- 1 / truncNuisanceWeight(NuisanceDenom = NuisanceDenom, 
                                                  MinNuisance = MinNuisance, 
                                                  RegimeName = names(PropScores)[a])
        return(list("PropScore" = PropScores[[a]],
                    "Hazards" = HazSurvPreds[[a]][["Hazards"]],
                    "EvntFreeSurv" = HazSurvPreds[[a]][["Survival"]][["TotalSurv"]],
                    "NuisanceWeight" = NuisanceWeight))
    })
    
    names(InitialEstimates) <- names(Regime)
    attr(InitialEstimates, "Times") <- Hazards[["Time"]]
    attr(InitialEstimates, "InitFits") <- InitFits
    return(InitialEstimates)
}

#' Truncate Nuisance Weights for Numerical Stability
#'
#' Applies lower bound truncation to the g-related nuisance weights (propensity
#' score times censoring survival) to prevent numerical instability from near-zero
#' denominators in the efficient influence curve.
#'
#' @param NuisanceDenom matrix (times x subjects); the product g(A|W) * S_C(t|A,W)
#' @param MinNuisance numeric in (0, 1); minimum allowed value. Values below this
#'   are set to MinNuisance. Typical default: 5 / (sqrt(n) * log(n))
#' @param RegimeName character; name of the intervention regime (for messages)
#'
#' @return matrix with same dimensions as NuisanceDenom, with values truncated
#'   at MinNuisance. Has attributes:
#'   \itemize{
#'     \item original: the pre-truncation values (if any truncation occurred)
#'     \item message: summary of truncation applied
#'   }
#'
#' @details
#' Truncation addresses the practical positivity assumption: even if theoretical
#' positivity holds, near-zero propensity scores or censoring survival can cause
#' the EIC to have very high variance. The bias-variance tradeoff controlled by
#' MinNuisance is fundamental to TMLE performance.
#'
#' A large percentage of truncated weights suggests potential positivity violations:
#' \itemize{
#'   \item >10% truncated: Consider whether the intervention is well-supported
#'   \item >25% truncated: Results may be unreliable; consider redefining the intervention
#'   \item >50% truncated: Severe positivity issues; results should not be trusted
#' }
#'
#' @keywords internal
truncNuisanceWeight <- function(NuisanceDenom, MinNuisance, RegimeName) {
    if (is.function(MinNuisance)) {
        warning("MinNuisance functions are not yet supported. Using default: ",
                "5 / (sqrt(n) * log(n)).", call. = FALSE)
        MinNuisance <- 5 / log(ncol(NuisanceDenom)) / sqrt(ncol(NuisanceDenom))
    }
    if (is.numeric(MinNuisance) & length(MinNuisance) == 1) {
        if (MinNuisance < 1 & MinNuisance > 0) {
            if (min(NuisanceDenom) < MinNuisance) {
                pct_subjects <- round(mean(apply(NuisanceDenom, 2, function(subj) any(subj < MinNuisance)), 3) * 100)
                pct_total <- round(mean(NuisanceDenom < MinNuisance), 3) * 100
                PositivityWarning <- paste0(
                    "For Intervention \"", RegimeName, "\": ",
                    pct_subjects, "% of subjects had at least one nuisance weight < ",
                    signif(MinNuisance, 3), "; ", pct_total, "% of all weights truncated")

                # Add severity guidance
                if (pct_subjects > 10) {
                    PositivityWarning <- paste0(PositivityWarning,
                        ". WARNING: High truncation rate suggests potential positivity violations.")
                }

                attr(NuisanceDenom, "original") <- NuisanceDenom
                attr(NuisanceDenom, "message") <- PositivityWarning
                NuisanceDenom[NuisanceDenom < MinNuisance] <- MinNuisance
            } else {
                attr(NuisanceDenom, "message") <- paste0(
                    "For Intervention \"", RegimeName, "\": ",
                    "no nuisance weights below threshold (", signif(MinNuisance, 3), ")")
            }

        } else {
            warning("MinNuisance must be in (0, 1). Current value: ", MinNuisance,
                    ". Nuisance weights will NOT be truncated, which may cause ",
                    "numerical instability if propensity scores approach 0 or 1.",
                    call. = FALSE)
            attr(NuisanceDenom, "message") <- paste0(
                "MinNuisance=", MinNuisance, " invalid (must be in (0,1)). ",
                "Weights not truncated - results may be unstable.")
        }
    }
    return(NuisanceDenom)
}

# screenCovRanger <- function(Data, j, nVar =  10, min.node.size = 3, mtry = floor(sqrt(ncol(Data))), 
#                             write.forest = FALSE, oob.error = FALSE, importance = "impurity", ...) 
# {
#     if (!requireNamespace("ranger", quietly = TRUE)) {
#         stop("screenCovRanger requires the 'ranger' package")
#     }
#     IDCol <- attr(Data, "ID")
#     TrtCol <- attr(Data, "Treatment")
#     TimeCol <- attr(Data, "EventTime")
#     TypeCol <- attr(Data, "EventType")
#     
#     SurvFormula <- paste0("Surv(time=", TimeCol, ", event=", TypeCol, "==", j, ")~.")
#     FitRanger <- ranger::ranger(formula = as.formula(SurvFormula), 
#                                 data = Data[, .SD, .SDcols = setdiff(colnames(Data), c(TrtCol, IDCol))], 
#                                 min.node.size = min.node.size, write.forest = write.forest, 
#                                 mtry = mtry, oob.error = oob.error, importance = importance)
#     
#     CovSelected <- rank(-FitRanger$variable.importance, na.last = NA, ties.method = "first") <= nVar
#     return(CovSelected)
# }

