
#' Title
#'
#' @param EventTime
#' @param EventType
#' @param Treatment
#' @param CovDataTable
#' @param ID
#' @param TargetTimes
#' @param TargetEvents
#' @param Models
#' @param CVArgs
#' @param NumUpdateSteps
#' @param OneStepEps
#' @param PropScoreCutoff
#' @param Verbose
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
doConCRTmle <- function(EventTime, EventType, Treatment, CovDataTable,
                        ID = NULL, TargetTimes = sort(unique(EventTime)),
                        TargetEvents = NULL, Models, CVArgs = NULL, NumUpdateSteps = 25,
                        OneStepEps = 0.1, PropScoreCutoff = 0.05, Verbose = FALSE, ...)
{
  # check & format parameters  ----------------------------------------------------------------
  args <- checkConCRTmleArgs(EventTime, EventType, Treatment, CovDataTable,
                             ID, TargetTimes, TargetEvents, Models, CVArgs, NumUpdateSteps,
                             OneStepEps, Verbose)
  Data <- args$Data
  Events <- sort(unique(Data$Event))
  Events <- Events[Events > 0]

  # helper functions  -------------------------------------------------------------------------

  PnEIC_norm_fun <- function(x, y) {
    return(sqrt(sum(unlist(x) * unlist(y))))
  }

  PnEIC_wt_fun <- function(PnEIC) {
    return(PnEIC)
  }

  # initial estimation ------------------------------------------------------------------------
  InitEsts <- getInitialEstimates(Data, CovdataTable, Models, PropScoreCutoff, TargetTimes)

  # initialize clever covariate table ---------------------------------------------------------
  ClevCovTbl <- getInitialCleverCov(InitEsts[["PropScore"]], InitEsts[["SupLrnFits"]],
                                    InitEsts[["Hazards"]], TargetEvents, TargetTimes, Events)

  ## initial estimator (g-computation) ----
  A_Fjk <- expand.grid(list("J" = TargetEvents, "K" = TargetTimes))
  A_Fjk <- c("A", apply(A_Fjk, 1, function(r) paste0("F.j", r["J"], ".t", r["K"])))
  Psi_init <- ClevCovTbl[, lapply(.SD, mean), .SDcols = A_Fjk, by = "A"][, -1]

  formatPsi <- function(Psi) {
    # Psi_init <- ClevCovTbl[, mget(c("A", grep("Psi.+", colnames(ClevCovTbl), value = T)))]
    # Psi_init <- Psi_init[, lapply(.SD, unique), by = "A"]
    Psi <- melt(Psi, id.vars = "A")
    Psi[, c("event", "time") := tstrsplit(variable, "\\.(j|t)", keep = 2:3, )]
    Psi <- dcast(Psi, A + time ~ event, value.var = "value")
    # Psi[, S := do.call(sum, mget(Events)), by = c("A", "time")]
    Psi <- Psi[, lapply(.SD, as.numeric)][order(time, A)]
  }


  # Update step ---------------------------------------------------------------------------------

  ## EIC ----
  IC <- getEIC(ClevCovTbl, TargetTimes, TargetEvents, Events)

  PnEIC <- IC$pnEIC

  PnEIC_wtd <- PnEIC_wt_fun(PnEIC)

  ClevCovTbl <- merge(ClevCovTbl, PnEIC_wtd, by = "A")
  PnEIC_norm <- PnEIC_norm_fun(PnEIC, PnEIC_wtd)

  PnEIC <- as.data.table(cbind("A" = 1:0,
                               rbind(PnEIC[grep("a1", names(PnEIC))],
                                     PnEIC[grep("a0", names(PnEIC))])))
  setnames(PnEIC, 2:ncol(PnEIC),
           do.call(paste0, expand.grid("PnEIC.j", TargetEvents, ".t", TargetTimes)))

  ## one-step tmle loop (one-step) ----

  for (step in 1:NumUpdateSteps) {
    if (Verbose)
      cat("starting step", step, "with update epsilon =", OneStepEps, "\n")
    PnEIC_prev <- PnEIC
    PnEIC_wtd_prev <- PnEIC_wtd
    PnEIC_norm_prev <- PnEIC_norm

    ## make backups ----
    ## save current cause-specific hazards and clever covs in case the update
    ## makes things worse

    OldCols <- grep("(h\\.j)|(Cox\\.j)", colnames(ClevCovTbl), value = TRUE)
    ClevCovTbl.old <- ClevCovTbl[, mget(OldCols)]

    ## 5.1 calculate update step direction -------------------------------------

    getUpdateStepDir <- function(l, j, k) {

    }
    ClevCovTbl[, (paste0("delta.j", Events, ".dx ")) := lapply(Events, function(l) 0)]

    for (l in Events) { # loop over event hazards
      ClevCovTbl[, (paste0("delta.j", l, ".dx")) := 0] # fluctuation of causes-specific hazard
      for (j in TargetEvents) { # loop over target Events
        for (k in TargetTimes) { # loop over target times
          ClevCovTbl[, (paste0("delta.j", l, ".dx")) := get(paste0("delta.j", l, ".dx")) +
               (Time <= k) * Ht.g *
               get(paste0("h.j", j, ".l", l, ".t", k)) *
               get(paste0("PnEIC.j", j, ".t", k)) / PnEIC_norm]
        }
      }
    }

    ## 5.3 update cause specific hazards ----------------------------------------
    ClevCovTbl[, S.t := 1]
    for (l in Events) { # loop over competing Events
      ClevCovTbl[, (paste0("Cox.j", l)) := get(paste0("Cox.j", l)) *
           exp(OneStepEps * get(paste0("delta.j", l, ".dx")))]
      #### truncation cox fit +- 500 why ? -----
      ClevCovTbl[, (paste0("Cox.j", l)) := pmax(-500, pmin(500, get(paste0("Cox.j", l))))]
      # Does changing the order of event updates matter? -------------------------------
      ClevCovTbl[, S.t := S.t * exp(-cumsum(get(paste0("BaseHaz.j", l)) *
                                      get(paste0("Cox.j", l)))),
         by = c("ID", "A")]
      # when does cox go to Infinity? -------
      if (any(is.infinite(ClevCovTbl[, get(paste0("Cox.j", l))])))
        stop(paste0("Cox.fit for j=", j, " went to infinity. wtf happened?\n"))
      if (min(ClevCovTbl[, S.t]) <= 0)
        stop(paste0("Survival at some time is leq 0. wtf happened?\n"))
    }
    ClevCovTbl[, "LagS.t" := c(1, S.t[-.N]), by = c("ID", "A")]

    ## update clever covariates ----
    for (j in TargetEvents) {
      ClevCovTbl[, (paste0("F.j", j, ".t")) := cumsum(`LagS.t` *
                                                get(paste0("BaseHaz.j", j)) *
                                                get(paste0("Cox.j", j))),
         by = c("ID", "A")]
      for (k in TargetTimes) {
        # ClevCovTbl[, (paste0("S.t", k)) := S.t[Time == k], by = c("ID", "A")]
        ClevCovTbl[, fjtau := get(paste0("F.j", j, ".t"))[Time == k], by = c("ID", "A")]
        for (l in Events) {
          ClevCovTbl[, h.jlk := 1] #### ASK HELENE: why is 1 the default? ---------
          #### event component ----
          ClevCovTbl[, h.jlk := (l == j) - (fjtau - get(paste0("F.j", j, ".t"))) / S.t]

          ClevCovTbl[round(S.t, 8) == 0, h.jlk := 0 - (l == j)] ##### ASK HELENE WHY -------

          ClevCovTbl[, (paste0("h.j", j, ".l", l, ".t", k)) := h.jlk]
        } #### mean treated/control event j risk at target time k ----
        ClevCovTbl[, (paste0("Psi.j", j, ".t", k)) := mean(fjtau), by = "A"]
        ClevCovTbl[, (paste0("F.j", j, ".t", k)) := fjtau]
      }
    }
    ClevCovTbl <- ClevCovTbl[, -c("h.jlk", "fjtau")]
    ## recalculate eic -----
    EicTbl <- data.table(ID = Data$ID)
    for (k in TargetTimes) {
      for (j in TargetEvents) {
        ClevCovTbl[, EIC := 0]
        for (l in Events) {
          ClevCovTbl[, EIC := EIC + get(paste0("h.j", j, ".l", l, ".t", k)) * Ht.g *
               (Time <= k) * (Time <= T.tilde) * (Trt == A) *
               ((Time == T.tilde) * (Event == l) -
                  get(paste0("Cox.j", l)) * get(paste0("BaseHaz.j", l)))]
        }
        EicTbl <- merge(EicTbl,
                        dcast(ClevCovTbl[, .(EIC = sum(EIC),
                                     Fjt = unique(get(paste0("F.j", j, ".t", k))),
                                     Psijt = unique(get(paste0("Psi.j", j, ".t", k)))),
                                 by = c("ID", "A")][, EIC := EIC + Fjt - Psijt],
                              ID ~ A, value.var = "EIC"),
                        by = "ID")
        setnames(EicTbl, "0", paste0("EIC.a0.j", j, ".t", k))
        setnames(EicTbl, "1", paste0("EIC.a1.j", j, ".t", k))
      }
    }

    ic <- list("ic" = EicTbl,
               "pnEIC" = colMeans(EicTbl[, -c("ID")]),
               "seEIC" = sqrt(diag(var(EicTbl[, -c("ID")]))))

    PnEIC <- ic$pnEIC
    PnEIC <- as.data.table(cbind("A" = 1:0,
                                 rbind(PnEIC[grep("a1", names(PnEIC))],
                                       PnEIC[grep("a0", names(PnEIC))])))
    eica <- do.call(paste0, expand.grid("PnEIC.j", TargetEvents, ".t", TargetTimes))
    setnames(PnEIC, 2:ncol(PnEIC), eica)

    PnEIC_wtd <- PnEIC_wt_fun(PnEIC)
    PnEIC_norm <- PnEIC_norm_fun(PnEIC, PnEIC_wtd)
    if (Verbose)
      cat("Step", step, "mean EIC norm = ", PnEIC_norm, "\n")

    if (PnEIC_norm_prev <= PnEIC_norm) {
      step <- step - 1
      warning(paste0("update overshot! one-step update epsilon of ",
                     OneStepEps, " will be halved\n"))
      OneStepEps <- OneStepEps / 2
      ## Revert to previous cshaz & clevcovs if update made things worse
      ClevCovTbl[, (OldCols) := ClevCovTbl.old[, mget(OldCols)]]

      PnEIC <- PnEIC_prev
      PnEIC_wtd <- PnEIC_wtd_prev
      PnEIC_norm <- PnEIC_norm_prev
    } else {
      # update PnEIC columns in ClevCovTbl with updated values
      ht_eic_cols <- grep("PnEIC", colnames(ClevCovTbl), value = T)
      ClevCovTbl <- ClevCovTbl[, !..ht_eic_cols]
      ClevCovTbl[PnEIC, on = .(A = A), (ht_eic_cols) := mget(sprintf("i.%s", ht_eic_cols))]

      summ_eic <- melt(ic$ic, id.vars = "ID")
      summ_eic[, c("dummy", "A", "J", "T") := tstrsplit(variable, "\\.(a|j|t)")]
      summ_eic <- summ_eic[, -c("variable", "dummy")]
      summ_eic[, A := paste0("a", get("A"))]
      summ_eic[, J := paste0("j", get("J"))]
      summ_eic[, "T" := paste0("t", get("T"))]
      summ_eic <- dcast(summ_eic, ID + A + `T` ~ J, value.var = "value")[
        , S := do.call(sum, mget(paste0("j", Events))),
        by = c("ID", "A", "T")]
      summ_eic <- dcast(summ_eic, ID ~ ...,
                        value.var = c(paste0("j", Events), "S"),
                        sep = ".")

      onestep_stop <- abs(colMeans(summ_eic[, -"ID"])) <=
        sqrt(colMeans(summ_eic[, -"ID"]^2)) / ( sqrt(nrow(Data)) * log(nrow(Data)))

      if (Verbose)
        cat(signif(abs(colMeans(summ_eic[, -"ID"])) /
                     (sqrt(colMeans(summ_eic[, -"ID"]^2)) /
                        ( sqrt(nrow(Data)) * log(nrow(Data)))), 2), "\n")


      if (all(onestep_stop) | step == NumUpdateSteps) {
        if (Verbose)
          ifelse(all(onestep_stop),
                 message("converged at step ", step),
                 warning("Warning: Algorithm did not converge by step ", step))

        ## format output ----
        Psi_tmle <- ClevCovTbl[, mget(c("A", grep("Psi.+", colnames(ClevCovTbl), value = T)))]
        Psi_tmle <- Psi_tmle[, lapply(.SD, unique), by = "A"]
        Psi_tmle <- melt(Psi_tmle, id.vars = "A")
        Psi_tmle[, c("dummy", "event", "time") := tstrsplit(variable, "\\.(j|t)")]
        Psi_tmle <- Psi_tmle[, -c("variable", "dummy")]
        Psi_tmle <- dcast(Psi_tmle, A + time ~ event, value.var = "value")
        Psi_tmle[, S := 1 - do.call(sum, mget(paste0(Events))), by = c("A", "time")]
        Psi_tmle <- Psi_tmle[, lapply(.SD, as.numeric)][order(time, A)]

        tmle_ic <- summ_eic[, -c("ID")]
        tmle_se <- sqrt(diag(var(tmle_ic)) / nrow(Data))

        break
      }
    }
  }

  # output --------------------------------------------------------------------------------------

  # g-comp (sl estimate)
  # unadjusted cox model
  # tmle & ic

  return(list(estimates = list("tmle" = Psi_tmle,
                               "g-comp" = Psi_init),
              ic = tmle_ic,
              se = tmle_se))

}

## check argument formats, types, and sizes ----
checkConCRTmleArgs <- function(EventTime, EventType, Treatment, CovDataTable,
                               ID = NULL, TargetTimes = sort(unique(EventTime)),
                               TargetEvents = NULL, Models, CVArgs = NULL, NumUpdateSteps = 25,
                               OneStepEps = 0.1, Verbose = FALSE, ...)
{
  if (!data.table::is.data.table(CovDataTable)) {
    CovDataTable <- data.table::as.data.table(CovDataTable)
    warning("CovDataTable must be a data.table. We have attempted to convert ",
            "the CovDataTable argument into an object of data.table class.\n")
  }

  if (any(c("ID",  "Event", "Trt", "t") %in% colnames(CovDataTable)))
    stop("'ID', 'Event', 'Trt', and 't' are reserved column",
         "names. Rename covariate column names to avoid name collisions.\n")


  if (!all(sapply(list(EventTime, EventType, Treatment, TargetTimes, TargetEvents),
                  function(vec) is.numeric(vec) | is.null(vec))))
    stop("EventTime, EventType, and Treatment ",
         "arguments must be numeric vectors\n")

  if (is.null(ID))
    ID = seq_along(EventTime)

  Data <- try(
    data.table::data.table("ID" = ID,
                           "Time" = EventTime,
                           "Event" = EventType,
                           "Trt" = Treatment,
                           CovDataTable)
  )
  if ("try-error" %in% class(Data)) {
    warning("Failed to create data datatable. ",
            "Check data inputs; see function help page\n")
    return(Data)
  }

  ReservedColumns <- c('Time', 'Event', 'Trt')
  if (!is.null(CovTrtTime)) {
    ReservedColumns <- c(ReservedColumns, "t")
    Data[, `t` := Time]
  }

  Events <- sort(unique(Data$Event))
  Events <- Events[Events > 0]

  if (is.null(TargetEvents))
    TargetEvents <- Events
  else if (!all(TargetEvents %in% Events))
    stop("")

  # target time(s)
  if (max(TargetTimes) >= max(Data["Event" != 0, "Time"])) {
    TargetTimes <- TargetTimes[TargetTimes < max(Data[["Time"]])]
    warning(paste0("No Observed events at max target time:",
                   " truncating target time(s) to be at or ",
                   "before the last observed event time"))
  }

  # To do:
  # check if covariates too highly correlated with ID, Time, Event, or Trt
  # target event(s)
  # binary treatment
  # stopping criteria
  # ...

  return(list(Data = Data))
}

getInitialEstimates <- function(Data, CovdataTable, Models, PropScoreCutoff,
                               TargetTimes) {
  ## cross validation setup ----
  # stratifying cv so that folds are balanced for treatment assignment & outcomes
  StrataIDs <- as.numeric(factor(paste0(Data[["Trt"]], ":", Data[["Event"]])))
  CVFolds <- origami::make_folds(n = Data, fold_fun = origami::folds_vfold,
                                 strata_ids = StrataIDs)

  ## PropScore score ----
  TrtTask <- sl3::make_sl3_Task(
    data = Data[, -c("Time", "Event", "ID")],
    covariates = colnames(CovDataTable),
    outcome = "Trt"
  )

  TrtSL <- sl3::Lrnr_sl$new(learners = Models[["A"]],
                            folds = CVFolds)
  TrtFit <- TrtSL$train(TrtTask)
  PropScore <- TrtFit$predict()

  if (min(PropScore) < PropScoreCutoff | max(PropScore) > 1 - PropScoreCutoff) {
    warning("practical positivity violation likely, truncating PropScore +-", PropScoreCutoff, "\n")
    PropScore[PropScore < PropScoreCutoff] <- PropScoreCutoff
    PropScore[PropScore > (1 - PropScoreCutoff)] <- 1 - PropScoreCutoff
  }

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
    Hazards <- rbind(data.table(time = 0, hazard = 0),
                     suppressWarnings(setDT(basehaz(ModelFit, centered = TRUE))))
    colnames(Hazards) <- c("Time", "SupLrnHaz")
    return(list("fit" = ModelFit, "BaseHaz" = Hazards))
  })

  ## output initial estimates ----
  HazTimes <- unique(c(TargetTimes, Data[["Time"]]))
  HazTimes <- HazTimes[HazTimes <= max(TargetTimes)]

  Hazards <- data.table("Time" = c(0, HazTimes))[order(Time)]
  Hazards <- rbind(cbind(Hazards, A = 0), cbind(Hazards, A = 1))

  ## baseline hazards ----
  for (j in grep("\\d+", names(Models), value = T)) {
    Hazards <- merge(Hazards, SupLrnFits[[j]][["BaseHaz"]], by = "Time", all.x = T)
    Hazards[, SupLrnHaz := zoo::na.locf(SupLrnHaz), by = "A"]
    if (j == "0") {
      Hazards[, SupLrnHaz := c(0, SupLrnHaz[-.N]), by = "A"]
      setnames(Hazards, "SupLrnHaz", "LagCumHaz.C.t")
    } else {
      # hazards[, paste0("cumBaseHaz.j", j) := haz]
      Hazards[, SupLrnHaz := c(0, diff(SupLrnHaz)), by = "A"]
      setnames(Hazards, "SupLrnHaz", paste0("BaseHaz.j", j))
    }
  }
  return(list("Hazards" = Hazards, "PropScore" = PropScore, "SupLrnFits" = SupLrnFits))
}

# initialize clever covariates
getInitialCleverCov <- function(PropScore, SupLrnFits, Hazards, TargetEvents, TargetTimes, Events) {
  ## clever covariate data.table ClevCovTbl
  ClevCovTbl <- rbind(cbind(Data, A = 1), cbind(Data, A = 0))
  setcolorder(ClevCovTbl, c("ID", "Trt", "Time", "Event"))

  ## clever covariate A component ----
  ClevCovTbl[A == 1, Ht.a := 1 / PropScore]
  ClevCovTbl[A == 0, Ht.a := 1 / (1 - PropScore)]

  ## clever covariate hazard components ----
  # A = target trt, Trt = assigned trt; Time = eval time, T.tilde = observed event time
  setnames(ClevCovTbl, c("A", "Trt"), c("Trt", "A"))
  for (j in grep("\\d+", names(Models), value = T)) {
    if (j == "0") {
      ClevCovTbl[, Cox.C := predict(SupLrnFits[["0"]]$fit, newdata = ClevCovTbl, type = "risk")]
    } else {
      ClevCovTbl[, (paste0("Cox.j", j)) := predict(SupLrnFits[[j]][["fit"]], newdata = ClevCovTbl, type = "risk")]
    }
  }
  setnames(ClevCovTbl, c("A", "Trt", "Time"), c("Trt", "A", "T.tilde"))

  ClevCovTblColNames <- c("ID", "Trt", "A", "T.tilde", "Event", "Ht.a",
                  grep("Cox", colnames(ClevCovTbl), value = T))
  ClevCovTbl <- ClevCovTbl[, ClevCovTblColNames, with = FALSE]
  ClevCovTbl <- merge(ClevCovTbl, Hazards, by = "A", allow.cartesian = T)[order(ID, Time, A)]
  setcolorder(ClevCovTbl, c("ID", "Time", "A"))

  ### clever covariate C component ----
  ClevCovTbl[, "LagS.C.a.t" := exp(-LagCumHaz.C.t * Cox.C)]
  if (min(ClevCovTbl[["LagS.C.a.t"]]) < 0.05) {
    warning(paste0("targeting a time where probability of being censored >95%,",
                   "cens survival truncated at min 0.05, but you should really",
                   "aim lower\n"))
  }
  ClevCovTbl[, "Ht.g" := Ht.a / LagS.C.a.t]
  ClevCovTbl <- ClevCovTbl[, -c("Ht.a", "Cox.C", "LagCumHaz.C.t", "LagS.C.a.t")]

  ### event-free survival ----
  ClevCovTbl[, S.t := 1]
  for (j in TargetEvents) {
    S.t_BaseHaz.j_Cox.j <- c("S.t", paste0("BaseHaz.j", j), paste0("Cox.j", j))
    ClevCovTbl[, S.t := .SD[["S.t"]] * exp(-cumsum(.SD[[2]] * .SD[[3]])),
       .SDcols = S.t_BaseHaz.j_Cox.j, by = c("ID", "A")]
    # ClevCovTbl[, S.t := S.t * exp(-cumsum(get(paste0("BaseHaz.j", j)) * get(paste0("Cox.j", j)))),
    #    by = c("ID", "A")]
  }
  if (min(ClevCovTbl[, S.t]) <= 0) {
    stop("you're targeting a time when there are no survivors, aim lower\n")
  }
  ClevCovTbl[, "LagS.t" := c(1, S.t[-.N]), by = c("ID", "A")]

  ### cause-specific risks ----
  for (j in TargetEvents) {
    LagS.t_BaseHaz.j_Cox.j <- c("LagS.t", paste0("BaseHaz.j", j), paste0("Cox.j", j))
    ClevCovTbl[, (paste0("F.j", j, ".t")) := cumsum(.SD[["LagS.t"]] * .SD[[2]] * .SD[[3]]),
       .SDcols = LagS.t_BaseHaz.j_Cox.j, by = c("ID", "A")]
  }

  ### cause-specific risks at target times ----
  ### target event (j), target time (k), cause (l) clever covariate ----
  for (k in TargetTimes) {
    for (j in TargetEvents) {
      #### event j risk at target time k ----
      Fjt_Time <- c(paste0("F.j", j, ".t"), "Time")
      ClevCovTbl[, (paste0("F.j", j, ".t", k)) := .SD[[1]][.SD[[2]] == k], .SDcols = Fjt_Time, by = c("ID", "A")]
      for (l in Events) {
        h.jlk <- paste0("h.j", j, ".l", l, ".t", k)
        Fjk_Fj <- c(paste0("F.j", j, ".t", k), paste0("F.j", j, ".t"))
        ClevCovTbl[, (h.jlk) := 1] # why was 1 the default, a 0 should squash oversights no?
        ClevCovTbl[, (h.jlk) := (l == j) - (.SD[[1]] - .SD[[2]]) / S.t, .SDcols = Fjk_Fj]
      }
    }
  }
  return(ClevCovTbl)
}

getEIC <- function(ClevCovTbl, TargetTimes, TargetEvents, Events) {
  EicTbl <- data.table(ID = unique(ClevCovTbl[["ID"]]))
  for (k in TargetTimes) {
    for (j in TargetEvents) {
      ClevCovTbl[, EIC := 0]
      for (l in Events) {
        EicCols <- c(paste0("h.j", j, ".l", l, ".t", k), paste0(c("Cox.j", "BaseHaz.j"), l),
                     "Event", "Time", "T.tilde", "Trt", "A", "Ht.g", "EIC")
        ClevCovTbl[(Time <= pmin(k, T.tilde)) & (Trt == A), EIC := EIC + Ht.g * .SD[[1]] *
                     ((Time == T.tilde) * (Event == l) - .SD[[2]] * .SD[[3]]),
                   .SDcols = EicCols]
      }
      EIC.jk <- ClevCovTbl[, list(EIC = sum(EIC), F.jk = mean(.SD[[2]])),
                           .SDcols = c("EIC", paste0("F.j", j, ".t", k)), by = c("ID", "A")]
      EIC.jk <- dcast(EIC.jk[, EIC := EIC + F.jk - mean(F.jk), by = "A"], ID ~ A, value.var = "EIC")
      setnames(EIC.jk, c("0", "1"), c(paste0("a0.j", j, ".t", k), paste0("a1.j", j, ".t", k)))

      EicTbl <- merge(EicTbl, EIC.jk, by = "ID")
    }
  }
  return(list("ic" = EicTbl[, -c("ID")],
              "pnEIC" = colMeans(EicTbl[, -c("ID")]),
              "seEIC" = sqrt(diag(var(EicTbl[, -c("ID")])))))
}
