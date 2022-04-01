#' Title
#'
#' @param Data data.table
#' @param Models list
#' @param CVFolds list
#' @param Hazards list
# #' @param HazFits list
# #' @param MinNuisance numeric (in the future a function)
# #' @param TargetEvents numeric vector
# #' @param TargetTimes numeric vector
# #' @param RegsOfInterest list
# #' @param Censored boolean
#'

getHazFit <- function(Data, Models, CVFolds, Hazards) {
    Time <- Event <- FitLP <- AtRisk <- basehaz <- BaseHaz <- NULL
    HazModels <- Models[grep("\\d+", names(Models))]
    HazModels <- lapply(seq_along(HazModels), function(j) {
        HazModel <- HazModels[[j]]
        attr(HazModel, "j") <- as.numeric(gsub(".*(\\d+).*", "\\1", names(HazModels[j])))
        return(HazModel)
    })
    SupLrnModels <- lapply(HazModels, function(Models_j) {
        SupLrnLibRisk <- data.table::data.table(matrix(NaN, nrow = nrow(Data), ncol = length(Models_j)))
        colnames(SupLrnLibRisk) <- names(Models_j)
        j <- attr(Models_j, "j")

        for (Fold_v in CVFolds) {
            TrainIndices <- Fold_v[["training_set"]]
            ValidIndices <- Fold_v[["validation_set"]]
            TrainData <- Data[TrainIndices, -c("ID")]
            ValidData <- Data[ValidIndices, ][order(-Time)]

            for (i in seq_along(Models_j)) {
                ## train model ----
                CoxphArgs <- list("formula" = Models_j[[i]], "data" = TrainData)
                ModelFit <- do.call(survival::coxph, CoxphArgs)

                ## validation loss (-log partial likelihood) ----
                ValidData[, FitLP := stats::predict(ModelFit, type = "lp", newdata = ValidData)]
                ValidData[, AtRisk := cumsum(exp(FitLP))]
                ValidData[AtRisk == 0, AtRisk := 1]
                ValidData[, names(Models_j)[i] := (Event == j) * (FitLP - log(AtRisk))]
            }
            SupLrnLibRisk[ValidIndices, names(Models_j) := subset(ValidData, select = names(Models_j))]
        }
        ## metalearner (discrete selector) ----
        SLCVRisk <- -colSums(SupLrnLibRisk)
        SLModel <- Models_j[[which.min(SLCVRisk)]]

        return(list("SupLrnCVRisks" = SLCVRisk, "SupLrnModel" = SLModel, "j" = j))
    })

    ## fit sl selection on full data ----
    HazFits <- lapply(SupLrnModels, function(SLMod) {
        ## fit sl model ----
        ModelFit <- do.call(survival::coxph, list("formula" = SLMod$SupLrnModel, "data" = Data))
        ## selected model must not contain interactions without including the separate terms also, e.g.
        #   we can't have ~ trt:sex without trt + sex. easiest way is to force e.g. trt*sex
        Hazards.j <- rbind(data.table(time = 0, hazard = 0),
                           suppressWarnings(setDT(basehaz(ModelFit, centered = TRUE))))
        colnames(Hazards.j) <- c("Time", "BaseHaz")
        Hazards.j <- merge(Hazards, Hazards.j, by = "Time", all.x = T)
        Hazards.j[, BaseHaz := c(0, diff(zoo::na.locf(BaseHaz)))]

        HazFitOut <- list("HazFit" = ModelFit, "BaseHaz" = Hazards.j, "HazModel" = SLMod)
        attr(HazFitOut, "j") <- SLMod[["j"]]
        return(HazFitOut)
    })
    names(HazFits) <- grep("\\d+", names(Models), value = TRUE)
    return(HazFits)
}

getHazSurvPred <- function(Data, HazFits, MinNuisance, TargetEvents,
                            TargetTimes, RegsOfInterest, Censored) {
    Targets <- expand.grid("Time" = TargetTimes, "Event" = TargetEvents)
    PredHazSurv <- lapply(RegsOfInterest, function(Reg) {
        PredData <- as.data.table(Data)[, "Trt" := Reg]
        PredHaz <- lapply(HazFits, function(HazFit) {
            exp.coef <- stats::predict(HazFit[["HazFit"]], newdata = PredData, type = "risk")
            haz <- sapply(exp.coef, function(expLP) HazFit[["BaseHaz"]][["BaseHaz"]]*expLP)
            attr(haz, "j") <- attr(HazFit, "j")
            return(haz)
        })
        names(PredHaz) <- names(HazFits)

        CensInd <- which(sapply(PredHaz, function(haz) attr(haz, "j") == 0))
        TotalSurv <- apply(do.call(`+`, PredHaz[-CensInd]), 2, function(haz) exp(-cumsum(haz)))

        if (Censored) {
            LaggedCensSurv <- apply(PredHaz[[CensInd]], 2, function(haz) c(1, utils::head(exp(-cumsum(haz)), -1)))
            PredHaz <- PredHaz[-CensInd]
        } else {
            LaggedCensSurv <- 1
        }

        Survival <- list("TotalSurv" = TotalSurv, "LaggedCensSurv" = LaggedCensSurv)
        return(list("Hazards" = PredHaz, "Survival" = Survival))
    })
    return(PredHazSurv)
}
