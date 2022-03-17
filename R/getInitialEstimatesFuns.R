getInitialEstimates <- function(Data, CovdataTable, Models, MinNuisanceDenom,
                                TargetTimes, RegsOfInterest, PropScoreBackend) {
    ## cross validation setup ----
    # stratifying cv so that folds are balanced for treatment assignment & outcomes
    # theory? but regressions may fail in practice with rare events otherwise
    StrataIDs <- as.numeric(factor(paste0(Data[["Trt"]], ":", Data[["Event"]])))
    CVFolds <- origami::make_folds(n = Data, fold_fun = origami::folds_vfold,
                                   strata_ids = StrataIDs)

    ## Propensity Scores for Regimes of Interest ----
    PropScores <- getPropScores(Data, CovDataTable, Models, MinNuisanceDenom, RegsOfInterest,
                                PropScoreBackend, CVFolds)

    ## hazards: Events & censoring ----
    SupLrnModels <- list()
    for (j in grep("\\d+", names(Models), value = T)) {
        Models_j <- Models[[j]]
        SupLrnLibRisk <- data.table::data.table(matrix(NaN, nrow = nrow(Data), ncol = length(Models_j)))
        colnames(SupLrnLibRisk) <- names(Models_j)

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
        SupLrnModels[[j]] <- list("SupLrnCVRisks" = -colSums(SupLrnLibRisk),
                                  "SupLrnModel" = Models_j[[which.min(-colSums(SupLrnLibRisk))]])
    }

    ## fit sl selection ----

    SupLrnFits <- lapply(SupLrnModels, function(SLMod) {
        ## fit sl model ----
        ModelFit <- do.call(survival::coxph, list("formula" = SLMod$SupLrnModel, "data" = Data))
        ## selected model must not contain interactions without including the separate terms also, e.g.
        #   we can't have ~ trt:sex without trt + sex. easiest way is to force e.g. trt*sex
        Hazards <- rbind(data.table(time = 0, hazard = 0),
                         suppressWarnings(setDT(basehaz(ModelFit, centered = TRUE))))
        colnames(Hazards) <- c("Time", "SupLrnHaz")
        return(list("fit" = ModelFit, "BaseHaz" = Hazards))
    })

    ## baseline hazards for obs times + target times ----
    HazTimes <- unique(c(TargetTimes, Data[["Time"]]))
    HazTimes <- HazTimes[HazTimes <= max(TargetTimes)]
    Hazards <- data.table("Time" = c(0, HazTimes))[order(Time)]

    ## baseline CumHaz for Censoring, haz for events ----
    for (j in grep("\\d+", names(Models), value = T)) {
        Hazards <- merge(Hazards, SupLrnFits[[j]][["BaseHaz"]], by = "Time", all.x = T)
        Hazards[, SupLrnHaz := zoo::na.locf(SupLrnHaz)]
        if (j == "0") {
            Hazards[, SupLrnHaz := c(0, SupLrnHaz[-.N])]
            setnames(Hazards, "SupLrnHaz", "LagCumHaz.C.t")
        } else {
            # hazards[, paste0("cumBaseHaz.j", j) := haz]
            Hazards[, SupLrnHaz := c(0, diff(SupLrnHaz))]
            setnames(Hazards, "SupLrnHaz", paste0("BaseHaz.j", j))
        }
    }
    return(list("Hazards" = Hazards[], "PropScores" = PropScores, "SupLrnFits" = SupLrnFits))
}
