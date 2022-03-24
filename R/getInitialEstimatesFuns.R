getInitialEstimates <- function(Data, CovDataTable, Models, MinNuisance, TargetEvents,
                                TargetTimes, RegsOfInterest, PropScoreBackend, Censored) {
    ## cross validation setup ----
    # stratifying cv so that folds are balanced for treatment assignment & outcomes
    # theory? but regressions may fail in practice with rare events otherwise
    StrataIDs <- as.numeric(factor(paste0(Data[["Trt"]], ":", Data[["Event"]])))
    CVFolds <- origami::make_folds(n = Data, fold_fun = origami::folds_vfold,
                                   strata_ids = StrataIDs)

    ## Propensity Scores for Regimes of Interest ----
    PropScores <- getPropScores(Data, CovDataTable, Models, MinNuisance, RegsOfInterest,
                                PropScoreBackend, CVFolds)
    TrtFit <- PropScores[["TrtFit"]]
    PropScores <- PropScores[["PropScores"]]

    ## hazards: Events & censoring ----
    ## baseline hazards for obs times + target times ----
    HazTimes <- unique(c(TargetTimes, Data[["Time"]]))
    HazTimes <- HazTimes[HazTimes <= max(TargetTimes)]
    Hazards <- data.table("Time" = c(0, HazTimes))[order(Time)]

    HazFits <- getHazFits(Data, Models, CVFolds, Hazards)
    HazSurvPreds <- getHazSurvPreds(Data, HazFits, MinNuisance, TargetEvents,
                                    TargetTimes, RegsOfInterest, Censored)
    InitialEstimates <- lapply(1:length(PropScores), function(a) {
        NuisanceWeight <- sapply(1:length(PropScores[[a]]), function(i) {
            PropScores[[a]][i] * HazSurvPreds[[a]][["Survival"]][["LaggedCensSurv"]][, i]})
        NuisanceWeight <- 1 / truncNuisanceDenom(NuisanceWeight, MinNuisance)
        return(list("PropScore" = PropScores[[a]],
                    "Hazards" = HazSurvPreds[[a]][["Hazards"]],
                    "EvntFreeSurv" = HazSurvPreds[[a]][["Survival"]][["TotalSurv"]],
                    "NuisanceWeight" = NuisanceWeight))
    })

    names(InitialEstimates) <- names(RegsOfInterest)
    attr(InitialEstimates, "times") <- Hazards[["Time"]]
    return(InitialEstimates)
}

getHazFits <- function(Data, Models, CVFolds, Hazards) {
    HazModels <- Models[grep("\\d+", names(Models))]
    HazModels <- lapply(1:length(HazModels), function(j) {
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

            for (i in 1:length(Models_j)) {
                ## train model ----
                CoxphArgs <- list("formula" = Models_j[[i]], "data" = TrainData)
                ModelFit <- do.call(survival::coxph, CoxphArgs)

                ## validation loss (-log partial likelihood) ----
                ValidData[, FitLP := predict(ModelFit, type = "lp", newdata = ValidData)]
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

getHazSurvPreds <- function(Data, HazFits, MinNuisance, TargetEvents,
                            TargetTimes, RegsOfInterest, Censored) {
    Targets <- expand.grid("Time" = TargetTimes, "Event" = TargetEvents)
    PredHazSurv <- lapply(RegsOfInterest, function(Reg) {
        PredData <- as.data.table(Data)[, "Trt" := Reg]
        PredHaz <- lapply(HazFits, function(HazFit) {
            exp.coef <- predict(HazFit[["HazFit"]], newdata = PredData, type = "risk")
            haz <- sapply(exp.coef, function(lambda) HazFit[["BaseHaz"]][["BaseHaz"]]*lambda)
            attr(haz, "j") <- attr(HazFit, "j")
            return(haz)
        })
        names(PredHaz) <- names(HazFits)

        CensInd <- which(sapply(PredHaz, function(haz) attr(haz, "j") == 0))
        TotalSurv <- apply(do.call(`+`, PredHaz[-CensInd]), 2, function(haz) exp(-cumsum(haz)))

        if (Censored) {
            LaggedCensSurv <- apply(PredHaz[[CensInd]], 2, function(haz) c(1, head(exp(-cumsum(haz)), -1)))
            PredHaz <- PredHaz[-CensInd]
        } else {
            LaggedCensSurv <- 1
        }

        Survival <- list("TotalSurv" = TotalSurv, "LaggedCensSurv" = LaggedCensSurv)
        return(list("Hazards" = PredHaz, "Survival" = Survival))
    })
    return(PredHazSurv)
}

