#' doConcrete
#'
#' @param ConcreteArgs "ConcreteArgs" sl3 object : output of formatArguments()
#'
# #' @param DataTable: data.table (N x ?)
# #' @param CovDataTable : data.table (N x ?)
# #' @param LongTime : numeric vector (?? x 1)
# #' @param ID : vector (N x 1)
# #' @param Events : numeric
# #' @param Censored : boolean
# #' @param TargetTime : numeric vector (length = K)
# #' @param TargetEvent : numeric vector \\subset EventType (length = J)
# #' @param Regime : list
# #' @param CVFolds : list
# #' @param Model : list of functions (length = L)
# #' @param PropScoreBackend : character
# #' @param HazEstBackend : character
# #' @param MaxUpdateIter : numeric
# #' @param OneStepEps : numeric
# #' @param MinNuisance : numeric
# #' @param Verbose : boolean
# #' @param GComp : boolean
# #' @param ReturnModels boolean
#'
#' @return object with s3 class "ConcreteEst"
#'
#' @importFrom grDevices dev.hold dev.flush devAskNewPage
#' @importFrom graphics abline
#' @importFrom stats density qnorm
#' @export doConcrete
#'
#' @examples
#' library(data.table)
#' library(survival)
#' data <- as.data.table(survival::pbc)
#' data[, trt := sample(0:1, nrow(data), TRUE)]
#' cols <- c("id", "time", "status", "trt",
#'           "age", "albumin", "sex", "bili")
#' data <- data[, .SD, .SDcols = cols]
#' 
#' intervention <- makeITT()
#' target.time <- 2500
#' target.event <- 1:2
#' model <- list("trt" = c("SL.glm", "SL.glmnet"),
#'               "0" = list(Surv(time, status == 0) ~ .),
#'               "1" = list(Surv(time, status == 1) ~ .),
#'               "2" = list(Surv(time, status == 2) ~ .))
#' 
#' # formatArguments() returns correctly formatted arguments for doConcrete()
#' concrete.args <- formatArguments(Data = data,
#'                                  EventTime = "time",
#'                                  EventType = "status",
#'                                  Treatment = "trt",
#'                                  ID = "id",
#'                                  Intervention = intervention,
#'                                  TargetTime = target.time,
#'                                  TargetEvent = target.event,
#'                                  Model = model,
#'                                  Verbose = FALSE)
#' 
#' # doConcrete() returns tmle (and g-comp plug-in) estimates of targeted risks
#' concrete.est <- doConcrete(concrete.args)
#' 
#' # getOutput returns risk difference, relative risk, and treatment-specific risks
#' concrete.out <- getOutput(Estimate = concrete.est, Estimand = c("rd", "rr", "risk"))
#' concrete.out$RD
#' concrete.out$RR
#' concrete.out$Risk

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

doConCRTmle <- function(DataTable, TargetTime, TargetEvent, Regime, CVFolds, Model, PropScoreBackend, 
                        HazEstBackend, MaxUpdateIter, OneStepEps, MinNuisance, Verbose, GComp, 
                        ReturnModels, ...)
{
    ratio <- Time <- Event <- PnEIC <- `seEIC/(sqrt(n)log(n))` <- NULL # for data.table compatibility w/ global var binding check
    
    # initial estimation ------------------------------------------------------------------------
    cat("Getting Initial Estimates:\n")
    Estimates <- getInitialEstimate(Data = DataTable, Model = Model, CVFolds = CVFolds, MinNuisance = MinNuisance,
                                    TargetEvent = TargetEvent, TargetTime = TargetTime, Regime = Regime,
                                    PropScoreBackend = PropScoreBackend, HazEstBackend = HazEstBackend, 
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
    cat("\nStarting TMLE Update:\n")
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
    class(Estimates) <- union("ConcreteEst", class(Estimates))
    print.ConcreteEst(Estimates, Verbose = Verbose)
    return(Estimates)
}

getNormPnEIC <- function(PnEIC, Sigma = NULL) {
    WeightedPnEIC <- PnEIC
    if (!is.null(Sigma)) {
        SigmaInv <- try(solve(Sigma))
        if (any(class(SigmaInv) == "try-error")) {
            SigmaInv <- solve(Sigma + diag(x = 1e-6, nrow = nrow(Sigma)))
            warning("regularization of Sigma needed for inversion")
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
    `Pt Est` <- se <- PnEIC <- `abs(PnEIC / Stop Criteria)` <- `seEIC/(sqrt(n)log(n))` <- NULL
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
        PnEICs[, "abs(PnEIC / Stop Criteria)" := abs(PnEIC / `seEIC/(sqrt(n)log(n))`)]
        PnEICs <- PnEICs[`abs(PnEIC / Stop Criteria)` > 1, .SD, 
                         .SDcols = c("Intervention", "Time", "Event", "Pt Est", "se", "PnEIC", 
                                     "abs(PnEIC / Stop Criteria)")]
        print(PnEICs[order(`abs(PnEIC / Stop Criteria)`, decreasing = TRUE), ], digits = 4)
        cat("\n")
    }
    
    for (a in seq_along(x)) {
        cat(attr(x[[a]]$NuisanceWeight, "message"), "\n")
    }
    cat("\n")
    
    cat("Initial Estimators:\n")
    TrtFit <- attr(x, "InitFits")[[1]]
    if (inherits(TrtFit, "SuperLearner")) {
        cat("Treatment: \n")
        if (is.matrix(TrtFit)) {
            print(TrtFit)
        } else {
            print(cbind(Risk = TrtFit$cvRisk, Coef = TrtFit$coef))
        }
    } else {
        cat("Treatment: \nPrinting for 'sl3' backend not yet enabled\n")
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
#' @param propscores logical: plot the densities of the propensity scores for each intervention
#' @param ask logical: whether or not to prompt for user input before displaying plots
#' @param ... additional arguments to be passed into plot methods
#' @exportS3Method plot ConcreteEst

plot.ConcreteEst <- function(x, convergence = TRUE, propscores = TRUE, ask = TRUE, ...) {
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if (convergence) {
        dev.hold()
        plot(x = seq_along(attr(x, "NormPnEICs")) - 1, y = attr(x, "NormPnEICs"), 
             xlab = "TMLE step", ylab = "PnEIC Norm", main = "TMLE Convergence", 
             pch = 20)
        dev.flush()
    }
    if (propscores) {
        for (i in seq_along(x)) {
            propscore <- x[[i]]$NuisanceWeight
            title = paste0("Distribution of Propensity Weights for Intervention \"",
                           names(x)[i], "\"")
            dev.hold()
            plot(density(1 / propscore), 
                 xlab = "G-related Nuisance Denominator", 
                 main = title)
            abline(v = 5 / sqrt(ncol(propscore)) / log(ncol(propscore)), col = "red")
            dev.flush()
        }
    }
}
