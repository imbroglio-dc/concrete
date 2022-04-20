#' Title
#'
#' @param Data data.table
#' @param Model list
#' @param CVFolds list
#' @param Hazards list
# #' @param HazFits list
# #' @param MinNuisance numeric (in the future a function)
# #' @param TargetEvent numeric vector
# #' @param TargetTime numeric vector
# #' @param Regime list
# #' @param Censored boolean
#'

getHazFit <- function(Data, EventTime, Model, CVFolds, Hazards) {
    Time <- Event <- FitLP <- AtRisk <- basehaz <- BaseHaz <- NULL
    HazModel <- Model[grep("\\d+", names(Model))]
    HazModel <- lapply(seq_along(HazModel), function(j) {
        haz.model <- HazModel[[j]]
        attr(haz.model, "j") <- as.numeric(gsub(".*(\\d+).*", "\\1", names(HazModel[j])))
        return(haz.model)
    })
    SupLrnModel <- lapply(HazModel, function(Model_j) {
        SupLrnLibRisk <- data.table::data.table(matrix(NaN, nrow = nrow(Data), ncol = length(Model_j)))
        colnames(SupLrnLibRisk) <- names(Model_j)
        j <- attr(Model_j, "j")
        id.col <- attr(Data, "ID")
        time.col <- attr(Data, "EventTime")
        event.col <- attr(Data, "EventType")

        for (Fold_v in CVFolds) {
            TrainIndices <- Fold_v[["training_set"]]
            ValidIndices <- Fold_v[["validation_set"]]
            TrainData <- Data[TrainIndices, .SD, .SDcols = setdiff(colnames(Data), id.col)]
            ValidData <- Data[ValidIndices, .SD, .SDcols = setdiff(colnames(Data), id.col)]
            setorderv(ValidData, cols = time.col, order = -1)

            for (i in seq_along(Model_j)) {
                ## train model ----
                CoxphArgs <- list("formula" = Model_j[[i]], "data" = TrainData)
                ModelFit <- do.call(survival::coxph, CoxphArgs)

                ## validation loss (-log partial likelihood) ----
                ValidData[, FitLP := stats::predict(ModelFit, type = "lp", newdata = ValidData)]
                ValidData[, AtRisk := cumsum(exp(FitLP))]
                ValidData[AtRisk == 0, AtRisk := 1]
                ValidData[, names(Model_j)[i] := (.SD == j) * (FitLP - log(AtRisk)), .SDcols = event.col]
            }
            SupLrnLibRisk[ValidIndices, names(Model_j) := subset(ValidData, select = names(Model_j))]
        }
        ## metalearner (discrete selector) ----
        SLCVRisk <- -colSums(SupLrnLibRisk)
        SLModel <- Model_j[[which.min(SLCVRisk)]]

        return(list("SupLrnCVRisks" = SLCVRisk, "SupLrnModel" = SLModel, "j" = j))
    })

    ## fit sl selection on full data ----
    HazFits <- lapply(SupLrnModel, function(SLMod) {
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
    names(HazFits) <- grep("\\d+", names(Model), value = TRUE)
    return(HazFits)
}

getHazSurvPred <- function(Data, HazFits, MinNuisance, TargetEvent,
                            TargetTime, Regime, Censored) {
    Target <- expand.grid("Time" = TargetTime, "Event" = TargetEvent)
    PredHazSurv <- lapply(Regime, function(Reg) {
        PredData <- as.data.table(Data)
        PredData[[attr(Data, "Treatment")]] <- Reg
        PredHaz <- lapply(HazFits, function(HazFit) {
            exp.coef <- stats::predict(HazFit[["HazFit"]], newdata = PredData, type = "risk")
            haz <- sapply(exp.coef, function(expLP) HazFit[["BaseHaz"]][["BaseHaz"]] * expLP)
            attr(haz, "j") <- attr(HazFit, "j")
            return(haz)(HazFit, "j")
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
