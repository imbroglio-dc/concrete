#' Title
#'
#' @param Data data.table
#' @param Model list
#' @param CVFolds list
#' @param Hazards list
#' @param ReturnModels boolean
# #' @param HazFits list
# #' @param MinNuisance numeric (in the future a function)
# #' @param TargetEvent numeric vector
# #' @param TargetTime numeric vector
# #' @param Regime list
# #' @param Censored boolean
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
        SupLrnLibRisk <- data.table::data.table(matrix(NaN, nrow = nrow(Data), ncol = length(ModelJ)))
        colnames(SupLrnLibRisk) <- names(ModelJ)
        j <- attr(ModelJ, "j")
        
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
        ## fit sl model ----
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
        
        HazFitOut <- list("HazFit" = ModelFit, "BaseHaz" = BaseHazJ)
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
    
    PredHazSurv <- lapply(Regime, function(Reg) {
        PredData <- as.data.table(Data)[, .SD, .SDcols = CovCols]
        setcolorder(PredData, neworder = CovCols)
        # PredData <- as.data.table(scale(PredData, center = TRUE, scale = FALSE))
        TrtNames <- colnames(Reg)
        PredData[, (TrtNames) := Reg[, .SD, .SDcols = TrtNames]]
        
        PredHaz <- lapply(HazFits, function(HazFit) {
            if (inherits(HazFit$HazFit, "cv.glmnet")) {
                PredData <- data.table::data.table(scale(PredData, center = TRUE, scale = FALSE))
                PredData[, (TrtNames) := Reg[, .SD, .SDcols = TrtNames]]
                exp.coef <- predict(HazFit$HazFit, newx = as.matrix(PredData), 
                                    s = HazFit$HazFit$lambda.min, type = "response")
                haz <- sapply(exp.coef, function(expLP) HazFit[["BaseHaz"]][["BaseHaz"]] * expLP)
            } else if (inherits(HazFit$HazFit, "coxph")) {
                exp.coef <- stats::predict(HazFit[["HazFit"]], newdata = PredData, type = "risk")
                haz <- sapply(exp.coef, function(expLP) HazFit[["BaseHaz"]][["BaseHaz"]] * expLP)
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

# SLCoxnet <- function(FitData, PredData, CovCols, TimeCol, TypeCol, j, alpha = 1, ...) {
#     if (!requireNamespace("glmnet", quietly = TRUE)) {
#         stop("SLCoxnet requires the 'glmnet' package")
#     }
#     FitLP <- AtRisk <- NULL
#     
#     ModelFit <- glmnet::glmnet(x = as.matrix(FitData)[, CovCols], 
#                                y = Surv(time = FitData[[TimeCol]], 
#                                         event = (FitData[[TypeCol]] == j), 
#                                         type = "right"),
#                                family = "cox", alpha = alpha, 
#                                penalty.factor = c(0, rep(1, length(CovCols) - 1)))
#     
#     for (s in seq_along(ModelFit$lambda)) {
#         PredTbl <- data.table()
#         PredTbl[, FitLP := stats::predict(object = ModelFit, s = ModelFit$lambda[s], type = "link", 
#                                           newx = as.matrix(PredData[, .SD, .SDcols = CovCols]))]
#         PredTbl[, AtRisk := cumsum(exp(FitLP))]
#         PredTbl[AtRisk == 0, AtRisk := 1]
#         PredTbl[, as.character(s) := (PredData[[TypeCol]] == j) * (FitLP - log(AtRisk))]
#     }
#     CoxnetRisk <- -colSums(PredTbl[, .SD, .SDcols = grep("coxnet", colnames(PredData), value = TRUE)])
#     LambdaMin <- ModelFit$lambda[as.numeric(names(which.min(CoxnetRisk)))]
#     
#     PredTbl[, names(ModelJ)[i] := get()]
#     PredTbl <- PredTbl[, .SD, .SDcols = !paste0("coxnet.s", seq_along(ModelFit$lambda))]
#     
#     return(list("CVRisk" = PredTbl, "ModelFit" = ModelFit))
# }
