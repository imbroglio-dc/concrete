#' Fit Cause-Specific Hazard Models
#'
#' Fits cause-specific hazard models for each event type using cross-validated
#' super learner ensemble methods. Supports Cox regression, elastic net (glmnet),
#' and survivalSL backends.
#'
#' @param Data data.table of observed data with attributes: EventTime, EventType,
#'   Treatment, ID, and CovNames specifying covariate column names
#' @param Model named list of model specifications, one for each unique EventType value.
#'   Each element can be: a formula (for Cox regression), a character vector of
#'   SuperLearner libraries, or a Lrnr.SurvivalSL specification from makeSurvivalSL()
#' @param CVFolds list of cross-validation fold assignments from origami::make_folds()
#' @param Hazards data.table with Time and Event columns specifying hazard evaluation times
#' @param ReturnModels logical; if TRUE, return fitted model objects for inspection
#'
#' @return list of class "HazFits" containing for each event type:
#'   - SupLrnCVRisks: cross-validated risks for each learner
#'   - SupLrnModel: the model specification used
#'   - j: the event type identifier
#'   - SLCoef: super learner coefficients
#'   - HazSL: fitted super learner model (if ReturnModels=TRUE)
#'
#' @keywords internal
#' @import survival
#' @importFrom stats predict
#'

getHazFit <- function(Data, Model, CVFolds, Hazards, ReturnModels) {
    Time <- Event <- FitLP <- AtRisk <- basehaz <- BaseHaz <- glmnet <- NULL

    IDCol <- attr(Data, "ID")
    TrtCol <- attr(Data, "Treatment")
    TimeCol <- attr(Data, "EventTime")
    TypeCol <- attr(Data, "EventType")
    CovDT <- subset(Data, select = attr(Data, "CovNames")[["ColName"]])
    HazModel <- Model[which(names(Model) %in% unique(Data[[TypeCol]]))]

    SupLrnModel <- lapply(HazModel, function(ModelJ) {
        j <- attr(ModelJ, "j")

        # Check if this is a survivalSL specification
        if (isSurvivalSLSpec(ModelJ)) {
            message("  Fitting survivalSL for event type ", j)
            SLFit <- fitSurvivalSL(
                Data = Data,
                j = j,
                TimeCol = TimeCol,
                TypeCol = TypeCol,
                TrtCol = TrtCol,
                IDCol = IDCol,
                methods = ModelJ$methods,
                metric = ModelJ$metric,
                cv = ModelJ$cv,
                seed = ModelJ$seed
            )

            if (is.null(SLFit)) {
                warning("survivalSL fit failed for event ", j,
                        ". Falling back to Cox model.")
                # Fall through to standard Cox fitting below
            } else {
                return(list(
                    "SupLrnCVRisks" = NULL,
                    "SupLrnModel" = ModelJ,
                    "j" = j,
                    "SLCoef" = SLFit$fit$weights$values,
                    "SurvivalSLFit" = SLFit
                ))
            }
        }
        # Standard Cox-based super learner fitting
        SupLrnLibRisk <- data.table::data.table(matrix(NaN, nrow = nrow(Data), ncol = length(ModelJ)))
        colnames(SupLrnLibRisk) <- names(ModelJ)

        for (Fold_v in CVFolds) {
            TrainIndices <- Fold_v[["training_set"]]
            ValidIndices <- Fold_v[["validation_set"]]
            TrainData <- Data[TrainIndices, .SD, .SDcols = setdiff(colnames(Data), IDCol)]
            ValidData <- Data[ValidIndices, .SD, .SDcols = setdiff(colnames(Data), IDCol)]
            setorderv(ValidData, cols = TimeCol, order = -1)
            
            ModelFits <- list()
            for (i in seq_along(ModelJ)) {
                if (!is.null(ModelJ[["screener"]])) {
                    if (ModelJ[["screener"]] == "ranger") {
                        message("screening not yet implemented")
                    } else 
                        if (is.function(ModelJ[["screener"]])) {
                            message("screening not yet implemented")
                        } else {
                            message("screening not yet implemented")
                        }
                }
                ## train model ----
                if (inherits(ModelJ[[i]], "Lrnr.Coxnet")) {
                    CovCols <- c(TrtCol, setdiff(colnames(TrainData), c(TimeCol, TypeCol, TrtCol, IDCol)))
                    ModelFit <- glmnet::glmnet(x = as.matrix(TrainData)[, CovCols], 
                                               y = Surv(time = TrainData[[TimeCol]], 
                                                        event = (TrainData[[TypeCol]] == j), 
                                                        type = "right"),
                                               family = "cox", 
                                               penalty.factor = c(0, rep(1, length(CovCols) - 1)))
                    z <- as.matrix(ValidData[, .SD, .SDcols = CovCols])
                } else { # Backend == "coxph"
                    CoxphArgs <- list("formula" = ModelJ[[i]], "data" = TrainData)
                    ModelFit <- do.call(survival::coxph, CoxphArgs)
                }
                if (ReturnModels) ModelFits[[i]] <- ModelFit
                
                ## validation loss (-log partial likelihood) ----
                if (inherits(ModelJ[[i]], "Lrnr.Coxnet")) {
                    for (s in seq_along(ModelFit$lambda)) {
                        ValidData[, FitLP := stats::predict(ModelFit, newx = z, s = ModelFit$lambda[s], type = "link")]
                        ValidData[, AtRisk := cumsum(exp(FitLP))]
                        ValidData[AtRisk == 0, AtRisk := 1]
                        ValidData[, paste0("coxnet.s", s) := (.SD == j) * (FitLP - log(AtRisk)), .SDcols = TypeCol]
                    }
                    CoxnetRisk <- -colSums(ValidData[, .SD, .SDcols = grep("coxnet", colnames(ValidData), value = TRUE)])
                    ValidData[, names(ModelJ)[i] := get(names(which.min(CoxnetRisk)))]
                    ValidData <- ValidData[, .SD, .SDcols = !paste0("coxnet.s", seq_along(ModelFit$lambda))]
                } else {
                    ValidData[, FitLP := stats::predict(ModelFit, type = "lp", newdata = ValidData)]
                    ValidData[, AtRisk := cumsum(exp(FitLP))]
                    ValidData[AtRisk == 0, AtRisk := 1]
                    ValidData[, names(ModelJ)[i] := (.SD == j) * (FitLP - log(AtRisk)), .SDcols = TypeCol]
                }
            }
            SupLrnLibRisk[ValidIndices, names(ModelJ) := subset(ValidData, select = names(ModelJ))]
        }
        ## metalearner (discrete selector) ----
        SLCVRisk <- -colSums(SupLrnLibRisk)
        SLModel <- ModelJ[[which.min(SLCVRisk)]]
        SLCoef <- as.numeric(SLCVRisk / min(SLCVRisk) == 1)
        names(SLCoef) <- names(SLCVRisk)
        
        if (ReturnModels) {
            return(list("SupLrnCVRisks" = SLCVRisk, "SupLrnModel" = SLModel, "j" = j, 
                        "SLCoef" = SLCoef, "ModelFits" = ModelFits))
        } else 
            return(list("SupLrnCVRisks" = SLCVRisk, "SupLrnModel" = SLModel, "j" = j, 
                        "SLCoef" = SLCoef))
    })
    names(SupLrnModel) <- sapply(SupLrnModel, function(sl) sl[["j"]])
    
    
    ## fit sl selection on full data ----
    HazFits <- lapply(SupLrnModel, function(SLMod) {
        ## Check if this is a survivalSL model
        if (!is.null(SLMod$SurvivalSLFit)) {
            # survivalSL model - already fitted during CV
            HazFitOut <- list(
                "HazFit" = SLMod$SurvivalSLFit,
                "BaseHaz" = NULL,  # Not used for survivalSL
                "isSurvivalSL" = TRUE
            )
            attr(HazFitOut, "j") <- SLMod[["j"]]
            attr(HazFitOut, "HazSL") <- list(
                "SupLrnCVRisks" = NULL,
                "j" = SLMod[["j"]],
                "SLCoef" = SLMod$SLCoef,
                "methods" = SLMod$SurvivalSLFit$methods
            )
            return(HazFitOut)
        }

        ## Standard Cox-based model fitting
        if (inherits(SLMod$SupLrnModel, "Lrnr.Coxnet")) {
            CovCols <- c(TrtCol, setdiff(colnames(Data), c(TimeCol, TypeCol, TrtCol, IDCol)))
            ModelFit <- glmnet::cv.glmnet(x = as.matrix(Data)[, CovCols],
                                          y = Surv(time = Data[[TimeCol]],
                                                   event = (Data[[TypeCol]] == SLMod[["j"]]),
                                                   type = "right"),
                                          family = "cox",
                                          penalty.factor = c(0, rep(1, length(CovCols) - 1)))
        } else {
            ## selected model must not contain interactions without including the separate terms also, e.g.
            #   we can't have ~ trt:sex without trt + sex. easiest way is to force e.g. trt*sex
            ModelFit <- do.call(survival::coxph, list("formula" = SLMod$SupLrnModel,
                                                      "data" = Data[, .SD, .SDcols = !attr(Data, "ID")]))
        }
        BaseHazCox <- paste0("Surv(time=", TimeCol, ", event=", TypeCol, "==", SLMod[["j"]], ")~",
                             paste0(TrtCol, collapse = "+"))
        BaseHazCox <- survival::coxph(as.formula(BaseHazCox),
                                      data = Data[, .SD, .SDcols = !attr(Data, "ID")])
        BaseHazJ <- rbind(data.table(time = 0, hazard = 0),
                          suppressWarnings(setDT(basehaz(BaseHazCox, centered = TRUE))))
        colnames(BaseHazJ) <- c("Time", "BaseHaz")
        BaseHazJ <- merge(Hazards, BaseHazJ, by = "Time", all.x = T)
        BaseHazJ[, BaseHaz := c(0, diff(zoo::na.locf(BaseHaz)))]

        HazFitOut <- list("HazFit" = ModelFit, "BaseHaz" = BaseHazJ, "isSurvivalSL" = FALSE)
        attr(HazFitOut, "j") <- SLMod[["j"]]
        SLMod[["SupLrnModel"]] <- NULL
        attr(HazFitOut, "HazSL") <- SLMod
        return(HazFitOut)
    })
    names(HazFits) <- names(SupLrnModel)
    return(HazFits)
}

getHazSurvPred <- function(Data, HazFits, MinNuisance, TargetEvent, TargetTime, Regime) {
    Censored <- any(Data[[attr(Data, "EventType")]] <= 0)
    Target <- expand.grid("Time" = TargetTime, "Event" = TargetEvent)
    IDCol <- attr(Data, "ID")
    TrtCol <- attr(Data, "Treatment")
    TimeCol <- attr(Data, "EventTime")
    TypeCol <- attr(Data, "EventType")
    CovCols <- c(TrtCol, setdiff(colnames(Data), c(TimeCol, TypeCol, TrtCol, IDCol)))

    # Get hazard times - must match the times grid from getInitialEstimate
    # This includes all observed event times + target times up to max target time
    HazTimes <- sort(unique(c(0, TargetTime, Data[[TimeCol]])))
    HazTimes <- HazTimes[HazTimes <= max(TargetTime)]

    # For Cox models, get times from BaseHaz if available
    if (!is.null(HazFits[[1]]$BaseHaz)) {
        HazTimes <- HazFits[[1]]$BaseHaz$Time
    }

    PredHazSurv <- lapply(Regime, function(Reg) {
        PredData <- as.data.table(Data)[, .SD, .SDcols = CovCols]
        setcolorder(PredData, neworder = CovCols)
        TrtNames <- colnames(Reg)
        PredData[, (TrtNames) := Reg[, .SD, .SDcols = TrtNames]]

        PredHaz <- lapply(HazFits, function(HazFit) {
            # Check if this is a survivalSL model
            if (isTRUE(HazFit$isSurvivalSL)) {
                SLFit <- HazFit$HazFit
                # Use the non-zero hazard times for survivalSL predictions
                SLTimes <- HazTimes[HazTimes > 0]
                haz <- predictSurvivalSLHazard(
                    SLFit = SLFit,
                    newdata = PredData,
                    times = SLTimes,
                    CovCols = SLFit$CovCols
                )
                # Prepend zero hazard row for time = 0
                haz <- rbind(rep(0, ncol(haz)), haz)
            } else if (inherits(HazFit$HazFit, "cv.glmnet")) {
                PredDataScaled <- data.table::data.table(scale(PredData, center = TRUE, scale = FALSE))
                PredDataScaled[, (TrtNames) := Reg[, .SD, .SDcols = TrtNames]]
                exp.coef <- predict(HazFit$HazFit, newx = as.matrix(PredDataScaled),
                                    s = HazFit$HazFit$lambda.min, type = "response")
                haz <- sapply(exp.coef, function(expLP) HazFit[["BaseHaz"]][["BaseHaz"]] * expLP)
            } else if (inherits(HazFit$HazFit, "coxph")) {
                exp.coef <- stats::predict(HazFit[["HazFit"]], newdata = PredData, type = "risk")
                haz <- sapply(exp.coef, function(expLP) HazFit[["BaseHaz"]][["BaseHaz"]] * expLP)
            } else {
                stop("Unknown HazFit type: ", class(HazFit$HazFit))
            }
            attr(haz, "j") <- attr(HazFit, "j")
            return(haz)
        })
        names(PredHaz) <- names(HazFits)
        
        CensInd <- which(sapply(PredHaz, function(haz) attr(haz, "j") <= 0))
        HazInd <- setdiff(seq_along(PredHaz), CensInd)

        TotalSurv <- apply(Reduce(`+`, PredHaz[HazInd]), 2, function(haz) exp(-cumsum(haz)))
        TotalSurv[TotalSurv < 1e-12] <- 1e-12

        if (Censored) {
            LaggedCensSurv <- apply(PredHaz[[CensInd]], 2, function(haz) c(1, utils::head(exp(-cumsum(haz)), -1)))
        } else {
            LaggedCensSurv <- 1
        }
        PredHaz <- PredHaz[HazInd]

        Survival <- list("TotalSurv" = TotalSurv, "LaggedCensSurv" = LaggedCensSurv)
        return(list("Hazards" = PredHaz, "Survival" = Survival))
    })
    return(PredHazSurv)
}

#' Fit a survivalSL model for cause-specific hazard estimation
#'
#' @param Data data.table with survival data
#' @param j numeric event type indicator
#' @param TimeCol character name of time column
#' @param TypeCol character name of event type column
#' @param TrtCol character name(s) of treatment column(s)
#' @param IDCol character name of ID column
#' @param methods character vector of survivalSL learner names
#' @param metric character loss function for survivalSL
#' @param cv numeric number of CV folds for survivalSL
#' @param seed numeric random seed
#' @param ... additional arguments passed to survivalSL
#'
#' @return list with fitted survivalSL model and metadata
#' @keywords internal
fitSurvivalSL <- function(Data, j, TimeCol, TypeCol, TrtCol, IDCol,
                          methods = c("LIB_AFTgamma", "LIB_PHgompertz"),
                          metric = "ibs", cv = 5, seed = NULL, ...) {
    if (!requireNamespace("survivalSL", quietly = TRUE)) {
        stop("fitSurvivalSL requires the 'survivalSL' package. ",
             "Install it with: install.packages('survivalSL')")
    }

    # Prepare data for survivalSL (needs data.frame, not data.table)
    CovCols <- setdiff(colnames(Data), c(TimeCol, TypeCol, IDCol))
    FitData <- as.data.frame(Data[, .SD, .SDcols = c(TimeCol, TypeCol, CovCols)])

    # Create binary event indicator for this event type
    FitData[[".event"]] <- as.numeric(FitData[[TypeCol]] == j)

    # Build formula
    FormulaTxt <- paste0("survival::Surv(", TimeCol, ", .event) ~ ",
                         paste0(CovCols, collapse = " + "))
    Formula <- stats::as.formula(FormulaTxt)

    # Fit survivalSL
    SLArgs <- list(
        formula = Formula,
        data = FitData,
        methods = methods,
        metric = metric,
        cv = cv,
        show_progress = FALSE
    )
    if (!is.null(seed)) SLArgs$seed <- seed

    SLFit <- tryCatch(
        do.call(survivalSL::survivalSL, SLArgs),
        error = function(e) {
            warning("survivalSL fitting failed for event ", j, ": ", e$message,
                    "\nFalling back to default Cox model.")
            return(NULL)
        }
    )

    if (is.null(SLFit)) return(NULL)

    return(list(
        "fit" = SLFit,
        "j" = j,
        "TimeCol" = TimeCol,
        "TypeCol" = TypeCol,
        "CovCols" = CovCols,
        "methods" = methods,
        "metric" = metric
    ))
}

#' Convert survivalSL survival predictions to hazard estimates
#'
#' @param SLFit fitted survivalSL object
#' @param newdata data.frame for prediction
#' @param times numeric vector of times for hazard estimation
#' @param CovCols character vector of covariate column names
#'
#' @return matrix of hazard estimates (times x subjects)
#' @keywords internal
predictSurvivalSLHazard <- function(SLFit, newdata, times, CovCols) {
    if (!requireNamespace("survivalSL", quietly = TRUE)) {
        stop("predictSurvivalSLHazard requires the 'survivalSL' package")
    }

    # Prepare prediction data
    PredData <- as.data.frame(newdata[, .SD, .SDcols = CovCols])

    # Get survival predictions at specified times
    Preds <- stats::predict(SLFit$fit, newdata = PredData, newtimes = times)

    # Extract super learner survival predictions (subjects x times)
    SurvMat <- Preds$predictions$sl

    # Convert survival to hazard: h(t) = -d(log S(t))/dt
    # For discrete: h(t_i) ≈ -log(S(t_i)/S(t_{i-1})) = log(S(t_{i-1})) - log(S(t_i))
    # Bound survival away from 0 to avoid log(0)
    SurvMat[SurvMat < 1e-10] <- 1e-10
    SurvMat[SurvMat > 1] <- 1

    LogSurv <- log(SurvMat)

    # Compute discrete hazard increments
    # Prepend 0 (log(1) = 0 for S(0) = 1)
    n_subj <- nrow(SurvMat)
    n_times <- ncol(SurvMat)

    HazMat <- matrix(0, nrow = n_times, ncol = n_subj)

    for (i in seq_len(n_subj)) {
        # First time point: h(t_1) = -log(S(t_1))
        HazMat[1, i] <- -LogSurv[i, 1]
        # Subsequent: h(t_k) = log(S(t_{k-1})) - log(S(t_k))
        if (n_times > 1) {
            HazMat[2:n_times, i] <- -diff(LogSurv[i, ])
        }
    }

    # Ensure non-negative hazards
    HazMat[HazMat < 0] <- 0

    return(HazMat)
}

#' Check if a model specification is for survivalSL
#'
#' @param ModelSpec model specification object
#' @return logical
#' @keywords internal
isSurvivalSLSpec <- function(ModelSpec) {
    inherits(ModelSpec, "Lrnr.SurvivalSL") ||
        (is.list(ModelSpec) && !is.null(ModelSpec$survivalSL))
}

#' Create a survivalSL model specification
#'
#' @param methods character vector of survivalSL learner names.
#'   Available learners include: "LIB_AFTgamma", "LIB_AFTggamma", "LIB_AFTllogis",
#'   "LIB_AFTweibull", "LIB_PHexponential", "LIB_PHgompertz", "LIB_PHspline",
#'   "LIB_RSF", "LIB_PLANN". Note: Cox-based learners (LIB_COX*) may have
#'   compatibility issues with some R versions.
#' @param metric character loss function for weight optimization.
#'   Options: "auc" (area under ROC), "ibs" (integrated Brier score),
#'   "ibll" (integrated binomial log-likelihood), "ll" (log-likelihood).
#' @param cv numeric number of cross-validation folds for survivalSL (default: 5)
#' @param seed numeric optional random seed for reproducibility
#' @param ... additional arguments passed to survivalSL::survivalSL
#'
#' @return a model specification object of class "Lrnr.SurvivalSL"
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a survivalSL specification with AFT and PH models
#' sl_spec <- makeSurvivalSL(
#'   methods = c("LIB_AFTgamma", "LIB_PHgompertz", "LIB_RSF"),
#'   metric = "ibs",
#'   cv = 5
#' )
#'
#' # Use in concrete model specification
#' model <- list(
#'   "trt" = c("SL.glm", "SL.glmnet"),
#'   "0" = sl_spec,  # censoring hazard
#'   "1" = sl_spec   # event hazard
#' )
#' }
makeSurvivalSL <- function(methods = c("LIB_AFTgamma", "LIB_PHgompertz"),
                           metric = "ibs",
                           cv = 5,
                           seed = NULL,
                           ...) {
    # Validate methods
    validMethods <- c("LIB_AFTgamma", "LIB_AFTggamma", "LIB_AFTllogis",
                      "LIB_AFTweibull", "LIB_COXaic", "LIB_COXall",
                      "LIB_COXen", "LIB_COXlasso", "LIB_COXridge",
                      "LIB_PHexponential", "LIB_PHgompertz", "LIB_PHspline",
                      "LIB_PLANN", "LIB_RSF")

    invalidMethods <- setdiff(methods, validMethods)
    if (length(invalidMethods) > 0) {
        warning("Unknown survivalSL methods will be ignored: ",
                paste(invalidMethods, collapse = ", "))
        methods <- intersect(methods, validMethods)
    }

    if (length(methods) < 2) {
        stop("survivalSL requires at least 2 learner methods")
    }

    # Validate metric
    validMetrics <- c("auc", "ibs", "ibll", "ribs", "ribll", "bll", "ll")
    if (!(metric %in% validMetrics)) {
        stop("Invalid metric. Choose from: ", paste(validMetrics, collapse = ", "))
    }

    spec <- list(
        survivalSL = TRUE,
        methods = methods,
        metric = metric,
        cv = cv,
        seed = seed,
        extra_args = list(...)
    )

    class(spec) <- c("Lrnr.SurvivalSL", "list")
    return(spec)
}
