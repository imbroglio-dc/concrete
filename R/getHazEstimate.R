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
#'

getHazFit <- function(Data, Model, CVFolds, Hazards, HazEstBackend, ReturnModels) {
    Time <- Event <- FitLP <- AtRisk <- basehaz <- BaseHaz <- NULL
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
            TimeCol <- attr(Data, "EventTime")
            TypeCol <- attr(Data, "EventType")
            
            for (Fold_v in CVFolds) {
                TrainIndices <- Fold_v[["training_set"]]
                ValidIndices <- Fold_v[["validation_set"]]
                TrainData <- Data[TrainIndices, .SD, .SDcols = setdiff(colnames(Data), IDCol)]
                ValidData <- Data[ValidIndices, .SD, .SDcols = setdiff(colnames(Data), IDCol)]
                setorderv(ValidData, cols = TimeCol, order = -1)
                
                ModelFits <- list()
                for (i in seq_along(ModelJ)) {
                    ## train model ----
                    CoxphArgs <- list("formula" = ModelJ[[i]], "data" = TrainData)
                    ModelFit <- do.call(survival::coxph, CoxphArgs)
                    if (ReturnModels) ModelFits[[i]] <- ModelFit
                    
                    ## validation loss (-log partial likelihood) ----
                    ValidData[, FitLP := stats::predict(ModelFit, type = "lp", newdata = ValidData)]
                    ValidData[, AtRisk := cumsum(exp(FitLP))]
                    ValidData[AtRisk == 0, AtRisk := 1]
                    ValidData[, names(ModelJ)[i] := (.SD == j) * (FitLP - log(AtRisk)), .SDcols = TypeCol]
                }
                SupLrnLibRisk[ValidIndices, names(ModelJ) := subset(ValidData, select = names(ModelJ))]
            }
            ## metalearner (discrete selector) ----
            SLCVRisk <- -colSums(SupLrnLibRisk)
            SLModel <- ModelJ[[which.min(SLCVRisk)]]
            
            if (ReturnModels) {
                return(c(list("SupLrnCVRisks" = SLCVRisk, "SupLrnModel" = SLModel, "j" = j), 
                         "ModelFits" = ModelFits))
            } else 
                return(list("SupLrnCVRisks" = SLCVRisk, "SupLrnModel" = SLModel, "j" = j))
        })
    } else {
        stop("Other hazard estimation methods not yet implemented")
    }
    
    
    ## fit sl selection on full data ----
    HazFits <- lapply(SupLrnModel, function(SLMod) {
        if (grepl("cox", tolower(HazEstBackend))) {
            ## fit sl model ----
            ModelFit <- do.call(survival::coxph, list("formula" = SLMod$SupLrnModel, "data" = Data))
            ## selected model must not contain interactions without including the separate terms also, e.g.
            #   we can't have ~ trt:sex without trt + sex. easiest way is to force e.g. trt*sex
            BaseHazJ <- rbind(data.table(time = 0, hazard = 0),
                              suppressWarnings(setDT(basehaz(ModelFit, centered = TRUE))))
            colnames(BaseHazJ) <- c("Time", "BaseHaz")
            BaseHazJ <- merge(Hazards, BaseHazJ, by = "Time", all.x = T)
            BaseHazJ[, BaseHaz := c(0, diff(zoo::na.locf(BaseHaz)))]
        }
        
        HazFitOut <- list("HazFit" = ModelFit, "BaseHaz" = BaseHazJ, "HazModel" = SLMod)
        attr(HazFitOut, "j") <- SLMod[["j"]]
        attr(HazFitOut, "HazSL") <- SLMod
        return(HazFitOut)
    })
    names(HazFits) <- grep("\\d+", names(Model), value = TRUE)
    return(HazFits)
}

getHazSurvPred <- function(Data, HazFits, MinNuisance, TargetEvent,
                           TargetTime, Regime, HazEstBackend) {
    Censored <- length(setdiff(TargetEvent, unique(Data[[attr(Data, "EventType")]]))) > 0
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
        
        CensInd <- which(sapply(PredHaz, function(haz) attr(haz, "j") == 0))
        HazInd <- setdiff(seq_along(PredHaz), CensInd)
        
        TotalSurv <- apply(Reduce(`+`, PredHaz[HazInd]), 2, function(haz) exp(-cumsum(haz)))
        
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
