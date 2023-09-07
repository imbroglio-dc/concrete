#' getInitialEstimate
#'
#' @param Data data.table
#' @param Model list
#' @param CVFolds : list
#' @param MinNuisance numeric
#' @param TargetEvent numeric vector
#' @param TargetTime numeric vector
#' @param Regime list
#' @param ReturnModels boolean

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

truncNuisanceWeight <- function(NuisanceDenom, MinNuisance, RegimeName) {
    if (is.function(MinNuisance)) {
        warning("MinNuisance functions are not yet supported.")
        MinNuisance <- 5 / log(ncol(NuisanceDenom)) / sqrt(ncol(NuisanceDenom))
    }
    if (is.numeric(MinNuisance) & length(MinNuisance) == 1) {
        if (MinNuisance < 1 & MinNuisance > 0) {
            if (min(NuisanceDenom) < MinNuisance) {
                PositivityWarning <- paste(
                    "For Intervention \"", RegimeName, "\", ",
                    round(mean(apply(NuisanceDenom, 2, function(subj) any(subj < MinNuisance)), 3) * 100),
                    "% of subjects had at least one G-related nuisance weight falling below ", 
                    signif(MinNuisance, 3), ", and ", 
                    round(mean(NuisanceDenom < MinNuisance), 3) * 100, "% of total G-related ", 
                    "nuisance weights were bounded to ", signif(MinNuisance, 3), sep = "")
                attr(NuisanceDenom, "original") <- NuisanceDenom
                attr(NuisanceDenom, "message") <- PositivityWarning
                NuisanceDenom[NuisanceDenom < MinNuisance] <- MinNuisance
            } else {
                attr(NuisanceDenom, "message") <- paste(
                    "For Intervention \"", RegimeName, "\", ", 
                    "no subjects had G-related nuisance weights falling below ", 
                    signif(MinNuisance, 3), sep = "")
            }
            
        } else {
            warning("MinNuisance improperly specified. G-related nuisance weights will not be ", 
                    "bounded away from 0, which can lead to computational instability.")
            attr(NuisanceDenom, "message") <- paste("MinNuisance improperly specified. G-related ",
                                                    "nuisance weights will not be bounded away from", 
                                                    " 0, which can lead to computational instability.", 
                                                    sep = "")
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

