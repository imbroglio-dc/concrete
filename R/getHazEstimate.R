#' Title
#'
#' @param Data data.table
#' @param Model list
#' @param CVFolds list
#' @param Hazards list
#' @param HazEstBackend : character
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

getHazFit <- function(Data, Model, CVFolds, Hazards, HazEstBackend, ReturnModels) {
    Time <- Event <- FitLP <- AtRisk <- basehaz <- BaseHaz <- glmnet <- NULL
    HazModel <- Model[grep("\\d+", names(Model))]
    HazModel <- lapply(seq_along(HazModel), function(j) {
        HazModelJ <- HazModel[[j]]
        attr(HazModelJ, "j") <- as.numeric(gsub(".*(\\d+).*", "\\1", names(HazModel[j])))
        return(HazModelJ)
    })
    if (grepl("cox", tolower(HazEstBackend))) {
        SupLrnModel <- lapply(HazModel, function(ModelJ) {
            SupLrnLibRisk <- data.table::data.table(matrix(NaN, nrow = nrow(Data), ncol = length(ModelJ)))
            colnames(SupLrnLibRisk) <- names(ModelJ)
            j <- attr(ModelJ, "j")
            IDCol <- attr(Data, "ID")
            TrtCol <- attr(Data, "Treatment")
            TimeCol <- attr(Data, "EventTime")
            TypeCol <- attr(Data, "EventType")
            CovDT <- subset(Data, select = attr(Data, "CovNames")[["ColName"]])
            
            for (Fold_v in CVFolds) {
                TrainIndices <- Fold_v[["training_set"]]
                ValidIndices <- Fold_v[["validation_set"]]
                TrainData <- Data[TrainIndices, .SD, .SDcols = setdiff(colnames(Data), IDCol)]
                ValidData <- Data[ValidIndices, .SD, .SDcols = setdiff(colnames(Data), IDCol)]
                setorderv(ValidData, cols = TimeCol, order = -1)
                
                ModelFits <- list()
                for (i in seq_along(ModelJ)) {
                    ## train model ----
                    if (ModelJ[[i]] == "coxnet") {
                        CovCols <- setdiff(colnames(ValidData), c(TimeCol, TypeCol, TrtCol, IDCol))
                        ModelFit <- glmnet(x = as.matrix(TrainData[, .SD, .SDcols = CovCols]), 
                                           y = Surv(time = TrainData[[TimeCol]], 
                                                    event = (TrainData[[TypeCol]] == j), 
                                                    type = "right"),
                                           family = "cox")
                        z <- as.matrix(ValidData[, .SD, .SDcols = CovCols])
                    } else {
                        CoxphArgs <- list("formula" = ModelJ[[i]], "data" = TrainData)
                        ModelFit <- do.call(survival::coxph, CoxphArgs)
                    }
                    if (ReturnModels) ModelFits[[i]] <- ModelFit
                    
                    ## validation loss (-log partial likelihood) ----
                    if (ModelJ[[i]] == "coxnet") {
                        for (s in seq_along(ModelFit$lambda)) {
                            ValidData[, FitLP := stats::predict(ModelFit, newx = as.matrix(z), s = ModelFit$lambda[s], type = "link")]
                            ValidData[, AtRisk := cumsum(exp(FitLP))]
                            ValidData[AtRisk == 0, AtRisk := 1]
                            ValidData[, paste0("coxnet.s", s) := (.SD == 1) * (FitLP - log(AtRisk)), .SDcols = "EVENT"]
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
    } else {
        stop("Other hazard estimation methods not yet implemented")
    }
    
    
    ## fit sl selection on full data ----
    HazFits <- lapply(SupLrnModel, function(SLMod) {
        if (grepl("cox", tolower(HazEstBackend))) {
            ## fit sl model ----
            if (SLMod$SupLrnModel == "coxnet") {
                
            } else {
                ## selected model must not contain interactions without including the separate terms also, e.g.
                #   we can't have ~ trt:sex without trt + sex. easiest way is to force e.g. trt*sex
                ModelFit <- do.call(survival::coxph, list("formula" = SLMod$SupLrnModel, 
                                                          "data" = Data[, .SD, .SDcols = !attr(Data, "ID")]))
                
                BaseHazJ <- rbind(data.table(time = 0, hazard = 0),
                                  suppressWarnings(setDT(basehaz(ModelFit, centered = TRUE))))
                colnames(BaseHazJ) <- c("Time", "BaseHaz")
                BaseHazJ <- merge(Hazards, BaseHazJ, by = "Time", all.x = T)
                BaseHazJ[, BaseHaz := c(0, diff(zoo::na.locf(BaseHaz)))]
            }
        }
        
        HazFitOut <- list("HazFit" = ModelFit, "BaseHaz" = BaseHazJ)
        attr(HazFitOut, "j") <- SLMod[["j"]]
        SLMod[["SupLrnModel"]] <- NULL
        attr(HazFitOut, "HazSL") <- SLMod
        return(HazFitOut)
    })
    names(HazFits) <- grep("\\d+", names(Model), value = TRUE)
    return(HazFits)
}

getHazSurvPred <- function(Data, HazFits, MinNuisance, TargetEvent,
                           TargetTime, Regime, HazEstBackend) {
    Censored <- any(Data[[attr(Data, "EventType")]] <= 0)
    Target <- expand.grid("Time" = TargetTime, "Event" = TargetEvent)
    PredHazSurv <- lapply(Regime, function(Reg) {
        PredData <- as.data.table(Data)
        PredData[[attr(Data, "Treatment")]] <- Reg
        PredHaz <- lapply(HazFits, function(HazFit) {
            exp.coef <- stats::predict(HazFit[["HazFit"]], newdata = PredData, type = "risk")
            haz <- sapply(exp.coef, function(expLP) HazFit[["BaseHaz"]][["BaseHaz"]] * expLP)
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
