#' Perform TMLE Estimation of Cause-Specific Absolute Risks
#'
#' Main estimation function for the concrete package. Computes continuous-time
#' one-step TMLE estimates of cause-specific cumulative incidence (absolute risk)
#' under specified intervention regimes.
#'
#' @param ConcreteArgs "ConcreteArgs" object; validated output from formatArguments().
#'   Contains all data, model specifications, and estimation parameters.
#'
#' @return Object of class "ConcreteEst" containing:
#'   \itemize{
#'     \item For each intervention regime: Hazards, EvntFreeSurv, PropScore,
#'       NuisanceWeight, SummEIC, and IC (efficient influence curve)
#'     \item Attributes: TmleConverged (list with convergence status and step count),
#'       NormPnEICs (trajectory of PnEIC norms), TargetTime, TargetEvent, T.tilde,
#'       Delta, GComp (whether g-computation estimates included)
#'   }
#'
#' @details
#' The estimation proceeds in three stages:
#' \enumerate{
#'   \item Initial estimation: Fits propensity scores and cause-specific hazards
#'     using cross-validated super learner ensembles
#'   \item EIC computation: Calculates the efficient influence curve for the
#'     cumulative incidence function under each intervention
#'   \item TMLE update: Iteratively updates hazard estimates along the least
#'     favorable submodel until the empirical mean of the EIC is sufficiently
#'     small (convergence criterion: |PnEIC| <= seEIC / (sqrt(n) * log(n)))
#' }
#'
#' The one-step TMLE update uses a small step size (controlled by OneStepEps)
#' and halves the step size if the update increases the PnEIC norm, ensuring
#' monotonic convergence.
#'
#' @importFrom grDevices dev.hold dev.flush devAskNewPage
#' @importFrom graphics abline
#' @importFrom stats density qnorm
#' @export doConcrete
#'
#' @examples
#' library(data.table)
#' library(concrete)
#' 
#' data <- as.data.table(survival::pbc)
#' data <- data[1:200, .SD, .SDcols = c("id", "time", "status", "trt", "age", "sex")]
#' data[, trt := sample(0:1, nrow(data), TRUE)]
#' 
#' # formatArguments() returns correctly formatted arguments for doConcrete()
#' 
#' concrete.args <- formatArguments(DataTable = data,
#'                                  EventTime = "time",
#'                                  EventType = "status",
#'                                  Treatment = "trt",
#'                                  ID = "id",
#'                                  TargetTime = 2500,
#'                                  TargetEvent = c(1, 2),
#'                                  Intervention = makeITT(),
#'                                  CVArg = list(V = 2))
#'                                  
#' # doConcrete() returns tmle (and g-formula plug-in) estimates of targeted risks
#' \donttest{concrete.est <- doConcrete(concrete.args)}

doConcrete <- function(ConcreteArgs) {
    if (!inherits(ConcreteArgs, "ConcreteArgs")) {
        stop("doConcrete takes a single argument, the output of the formatArguments() function. ", 
             "Run that function first and pass the resulting output into doConcrete() after ", 
             "addressing any errors or warnings.")
    }
    ArgList <- lapply(ls(ConcreteArgs), function(x) ConcreteArgs[[x]])
    names(ArgList) <- ls(ConcreteArgs)
    return(do.call(doConCRTmle, ArgList))
}

#' Internal TMLE Implementation
#'
#' Core implementation of the continuous-time TMLE algorithm. Called by doConcrete()
#' after argument validation.
#'
#' @param DataTable data.table; formatted data with attributes for column names
#' @param TargetTime numeric vector; time points at which to estimate risks
#' @param TargetEvent numeric vector; event types to target
#' @param Regime list; intervention specifications with g.star functions
#' @param CVFolds list; cross-validation fold assignments from origami
#' @param Model list; model specifications for propensity and hazard estimation
#' @param MaxUpdateIter integer; maximum TMLE update iterations
#' @param OneStepEps numeric; initial step size for TMLE updates (0, 1]
#' @param MinNuisance numeric; lower bound for propensity score truncation
#' @param Verbose logical; whether to print progress information
#' @param GComp logical; whether to compute g-computation estimates
#' @param ReturnModels logical; whether to return fitted model objects
#' @param ... additional arguments (unused, for compatibility)
#'
#' @return ConcreteEst object (see doConcrete documentation)
#' @keywords internal
doConCRTmle <- function(DataTable, TargetTime, TargetEvent, Regime, CVFolds, Model, MaxUpdateIter,
                        OneStepEps, MinNuisance, Verbose, GComp, ReturnModels, ...)
{
    ratio <- Time <- Event <- PnEIC <- `seEIC/(sqrt(n)log(n))` <- NULL # for data.table compatibility w/ global var binding check
    
    # initial estimation ------------------------------------------------------------------------
    Estimates <- getInitialEstimate(Data = DataTable, Model = Model, CVFolds = CVFolds,
                                    MinNuisance = MinNuisance, TargetEvent = TargetEvent,
                                    TargetTime = TargetTime, Regime = Regime,
                                    ReturnModels = ReturnModels)

    # get initial EIC (possibly with GComp plug-in estimate) ---------------------------------------------
    Estimates <- getEIC(Estimates = Estimates, Data = DataTable, Regime = Regime,
                        TargetEvent = TargetEvent, TargetTime = TargetTime,
                        MinNuisance = MinNuisance, GComp = GComp)
    
    
    # Update step -------------------------------------------------------------------------------
    ## Check if EIC is solved sufficienty and return outputs ----
    ## check PnEIC <= seEIC / (sqrt(n) log(n))
    SummEIC <- do.call(rbind, lapply(seq_along(Estimates), function(a) {
        cbind("Trt" = names(Estimates)[a], Estimates[[a]][["SummEIC"]])}))
    NormPnEIC <- getNormPnEIC(SummEIC[Time %in% TargetTime & Event %in% TargetEvent, PnEIC])
    OneStepStop <- SummEIC[, list("check" = abs(PnEIC) <= `seEIC/(sqrt(n)log(n))`,
                                  "ratio" = abs(PnEIC) / `seEIC/(sqrt(n)log(n))`),
                           by = c("Trt", "Time", "Event")]
    
    if (Verbose) printOneStepDiagnostics(OneStepStop, NormPnEIC)
    
    ## one-step tmle loop (one-step) ----
    message("\nStarting TMLE Update:\n")
    if (!all(sapply(OneStepStop[["check"]], isTRUE))) {
        Estimates <- doTmleUpdate(Estimates = Estimates, SummEIC = SummEIC, Data = DataTable,
                                  TargetEvent = TargetEvent, TargetTime = TargetTime,
                                  MaxUpdateIter = MaxUpdateIter, OneStepEps = OneStepEps,
                                  NormPnEIC = NormPnEIC, Verbose = Verbose)
    } else {
        attr(Estimates, "TmleConverged") <- list("converged" = TRUE, "step" = 0)
        attr(Estimates, "NormPnEICs") <- NormPnEIC
    }
    
    attr(Estimates, "TargetTime") <- TargetTime
    attr(Estimates, "T.tilde") <- DataTable[[attr(DataTable, "EventTime")]]
    attr(Estimates, "TargetEvent") <- TargetEvent
    attr(Estimates, "Delta") <- DataTable[[attr(DataTable, "EventType")]]
    attr(Estimates, "GComp") <- GComp

    # Note if survival warnings were generated (already attached from getInitialEstimate)
    if (!is.null(attr(Estimates, "SurvivalWarnings"))) {
        message("Note: Low survival warnings generated. Details: attr(result, 'SurvivalWarnings')")
    }

    class(Estimates) <- union("ConcreteEst", class(Estimates))
    return(Estimates)
}

#' Compute Norm of Empirical Mean EIC
#'
#' Calculates the (optionally weighted) Euclidean norm of the empirical mean
#' efficient influence curve. Used as the convergence criterion for TMLE.
#'
#' @param PnEIC numeric vector; empirical mean of the EIC for each target
#' @param Sigma matrix (optional); covariance matrix for weighted norm.
#'   If NULL, computes unweighted Euclidean norm.
#'
#' @return numeric; the norm ||PnEIC|| (or ||PnEIC||_Sigma if Sigma provided)
#'
#' @details
#' The weighted norm is computed as sqrt(PnEIC' * Sigma^{-1} * PnEIC).
#' If Sigma is singular, a small ridge regularization (1e-6) is added
#' to the diagonal before inversion.
#'
#' @keywords internal
getNormPnEIC <- function(PnEIC, Sigma = NULL) {
    WeightedPnEIC <- PnEIC
    if (!is.null(Sigma)) {
        SigmaInv <- try(solve(Sigma))
        if (any(class(SigmaInv) == "try-error")) {
            SigmaInv <- solve(Sigma + diag(x = 1e-6, nrow = nrow(Sigma)))
            warning("Covariance matrix Sigma was singular; added ridge regularization ",
                    "(1e-6) for inversion. This may occur with highly correlated EICs.")
        }
        WeightedPnEIC <- PnEIC %*% SigmaInv
    }
    return(sqrt(sum(unlist(PnEIC) * unlist(WeightedPnEIC))))
}

#' @describeIn doConcrete print.ConcreteEst print method for "ConcreteEst" class
#' @param x a ConcreteEst object
#' @param ... additional arguments to be passed into print methods
#' @exportS3Method print ConcreteEst
print.ConcreteEst <- function(x, ...) {
    `.` <- `..a` <- `Pt Est` <- se <- PnEIC <- `|Pn EIC| / Stop Criteria` <- `seEIC/(sqrt(n)log(n))` <- NULL
    cat("Continuous-Time One-Step TMLE targeting the Cause-Specific Absolute Risks for:\n")
    cat("Intervention", ifelse(length(x) > 1, "s", ""), ": ", 
        paste0("\"", names(x), "\"", collapse = ", "), "  |  ", sep = "")
    cat("Target Event", ifelse(length(attr(x, "TargetEvent")) > 1, "s", ""), ": ", 
        paste0(attr(x, "TargetEvent"), collapse = ", "), "  |  ", sep = "")
    cat("Target Time", ifelse(length(attr(x, "TargetTime")) > 1, "s", ""), ": ", 
        ifelse(length(attr(x, "TargetTime")) > 6,
               paste0(paste0(head(attr(x, "TargetTime"), 3), collapse = ", "), ",...,", 
                      paste0(tail(attr(x, "TargetTime"), 3), collapse = ", "), collapse = ""), 
               paste0(attr(x, "TargetTime"), collapse = ", ")), "\n\n", sep = "")
    
    cat(ifelse(isTRUE(attr(x, "TmleConverged")$converged), 
               paste0("TMLE converged at step ", attr(x, "TmleConverged")$step), 
               paste0("**TMLE did not converge!!**")), "\n\n")
    if (!(isTRUE(attr(x, "TmleConverged")$converged))) {
        PnEICs <- do.call(rbind, lapply(seq_along(x), function(a) 
            cbind("Intervention" = names(x)[a], x[[a]]$SummEIC)))
        Risks <- getRisk(x, TargetTime = attr(x, "TargetTime"), TargetEvent = attr(x, "TargetEvent"), 
                         GComp = FALSE)[, .SD, .SDcols = c("Intervention", "Time", "Event", "Pt Est", "se")]
        Risks <- rbind(Risks, Risks[, list("Event" = -1, "Pt Est" = 1 - sum(`Pt Est`), 
                                           "se" = sqrt(sum(se^2))), by = c("Intervention", "Time")])
        PnEICs <- merge(PnEICs, Risks, by = c("Intervention", "Time", "Event"))
        PnEICs[, "|Pn EIC| / Stop Criteria" := abs(PnEIC / `seEIC/(sqrt(n)log(n))`)]
        PnEICs <- PnEICs[`|Pn EIC| / Stop Criteria` > 1, .SD, 
                         .SDcols = c("Intervention", "Time", "Event", "Pt Est", "se", "PnEIC", 
                                     "|Pn EIC| / Stop Criteria")]
        print(PnEICs[order(`|Pn EIC| / Stop Criteria`, decreasing = TRUE), ], digits = 4)
        cat("\n")
    }
    
    for (a in seq_along(x)) {
        cat(attr(x[[a]]$NuisanceWeight, "message"), "\n")
    }
    cat("\n")
    
    cat("Initial Estimators:\n")
    for (a in setdiff(names(attr(x, "InitFits")), unique(attr(x, "Delta")))) {
        cat("Treatment \"", a, "\" :\n", sep = "")
        if (inherits(attr(x, "InitFits")[[a]], "SuperLearner")) {
            if (is.matrix(attr(x, "InitFits")[[a]])) {
                names(attr(x, "InitFits")[[a]]) <- sub(pattern = "Coef", replacement = "SL Weight", 
                                                       x = names(attr(x, "InitFits")[[a]]))
                print(attr(x, "InitFits")[[a]])
            } else {
                print(cbind(Risk = attr(x, "InitFits")[[a]]$cvRisk, 
                            "SL Weight" = attr(x, "InitFits")[[a]]$coef))
            }
            cat("\n")
        } else {
            cat("Treatment \"", a, "\": Printing for non-'SuperLearner' learners not yet enabled\n")
        }
    }
    
    cat("\n")
    for (Delta in sort(unique(attr(x, "Delta")))) {
        JFit <- attr(x, "InitFits")[[as.character(Delta)]]
        cat(ifelse(Delta <= 0, "Cens. ", "Event "), Delta, ": \n", sep = "")
        print(cbind(Risk = JFit$SupLrnCVRisks, Coef = JFit$SLCoef))
        cat("\n")
    }
}

#' @describeIn doConcrete plot.ConcreteEst plot method for "ConcreteEst" class
#' @param x a ConcreteEst object
#' @param convergence logical: plot the PnEIC norms for each TMLE small update step
#' @param gweights logical: plot the densities of the intervention-related nuisance weights for each intervention
#' @param ask logical: whether or not to prompt for user input before displaying plots
#' @param ... additional arguments to be passed into plot methods
#' @exportS3Method plot ConcreteEst

plot.ConcreteEst <- function(x, convergence = FALSE, gweights = TRUE, ask = FALSE, ...) {
    Intervention <- gDenomWeight <- y <- NULL
    if(!requireNamespace("ggplot2", quietly = TRUE))
        stop("Plotting requires the 'ggplot2' package")
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    fig <- list()
    if (convergence) {
        dev.hold()
        fig.conv <- ggplot2::ggplot(data.frame(x = seq_along(attr(x, "NormPnEICs")) - 1, 
                                               y = attr(x, "NormPnEICs"))) +
            ggplot2::geom_point(ggplot2::aes(x = x, y = y), size = 0) +
            ggplot2::labs(title = "TMLE Convergence", x = "TMLE step", y = "PnEIC Norm") + 
            ggplot2::theme_minimal()
        plot(fig.conv)
        fig <- c(fig, list("TMLEConvergence" = fig.conv))
        dev.flush()
    }
    if (gweights) {
        dev.hold()
        n <- ifelse(is.matrix(x[[1]]$NuisanceWeight), 
                    ncol(x[[1]]$NuisanceWeight), 
                    length(x[[1]]$NuisanceWeight))
        ps <- do.call(rbind, 
                      lapply(seq_along(x), function(a) {
                          g <- 1 / x[[a]]$NuisanceWeight
                          if (!is.null(attr(g, "original")))
                              g <- attr(g, "original")
                          data.frame("Intervention" = names(x)[a], 
                                     "gDenomWeight" = as.numeric(g))
                      }))
        fig.ps <- ggplot2::ggplot(ps) + ggplot2::lims(x = c(0, NA)) + ggplot2::theme_minimal() + 
            ggplot2::geom_line(mapping = ggplot2::aes(x = gDenomWeight, colour = Intervention), 
                               stat = "density", trim = TRUE) + 
            ggplot2::geom_vline(ggplot2::aes(xintercept = 5/(sqrt(n)*log(n))), colour = "red") +
            ggplot2::labs(title = "Distribution of Intervention-Related Nuisance Weights", 
                          subtitle = "Weights close to 0 warn of possible positivity violations", 
                          x = expression(~pi*"(a|w) "*S[c]*"(t|a,w)"), y = "Density")
        plot(fig.ps)
        fig <- c(fig, list("PropScores" = fig.ps))
        dev.flush()
    }
    invisible(fig)
}
