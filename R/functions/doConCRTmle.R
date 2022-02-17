
doConCRTmle <- function(EventTime, EventType, Treatment, CovDataTable, CovTrtTime = NULL,
                        ID = NULL, TargetTimes = sort(unique(EventTime)),
                        TargetEvents = NULL, Models, CVArgs = NULL, NumUpdateSteps = 25,
                        OneStepEps = 0.1, PropScoreCutoff = 0.05, Verbose = FALSE, ...)
{
  # parameter checking --------------------------------------------------------------------------

  ## check argument formats, types, and sizes ----
  checkConCRTmleArgs <- function(EventTime, EventType, Treatment, CovDataTable, CovTrtTime = NULL,
                                 ID = NULL, TargetTimes = sort(unique(EventTime)),
                                 TargetEvents = NULL, Models, CVArgs = NULL, NumUpdateSteps = 25,
                                 OneStepEps = 0.1, Verbose = FALSE, ...)
  {
    if (!data.table::is.data.table(CovDataTable)) {
      CovDataTable <- data.table::as.data.table(CovDataTable)
      warning("CovDataTable must be a data.table. We have attempted to convert ",
              "the CovDataTable argument into an object of data.table class.\n")
    }

    if (any(c("ID", "CovTrtTime", "Event", "Trt", "t") %in% colnames(CovDataTable)))
      stop("'ID', 'CovTrtTime', 'Event', 'Trt', and 't' are reserved column",
           "names. Rename covariate column names to avoid name collisions.\n")


    if (!all(sapply(list(EventTime, EventType, Treatment, CovTrtTime, TargetTimes, TargetEvents),
                    function(vec) is.numeric(vec) | is.null(vec))))
      stop("EventTime, EventType, Treatment, and CovTrtTime ",
           "arguments must be numeric vectors\n")

    ## check if covariates too highly correlated with ID, Time, CovTrtTime, Event, or Trt


    if (!is.null(CovTrtTime)) {
      if (is.null(ID))
        stop("When providing time-varying covariates and/or ",
             "treatments, an ID vector must be provided.\n")
      if (length(CovTrtTime) != nrow(CovDataTable))
        stop("CovTrtTime must be a numeric vector with length ",
             "equal to the number of rows in CovDataTable. \n")
      CovDataTable[, 't' := CovTrtTime]
    } else if (is.null(ID))
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

    # target event(s)
    # binary treatment
    # stopping criteria
    # ...

    return(list(Data = Data))
  }

  args <- checkConCRTmleArgs(EventTime, EventType, Treatment, CovDataTable, CovTrtTime,
                             ID, TargetTimes, TargetEvents, Models, CVArgs, NumUpdateSteps,
                             OneStepEps, Verbose)
  Data <- args$Data
  Events <- sort(unique(Data$Event))
  Events <- Events[Events > 0]


  # helper functions  ---------------------------------------------------------------------

  PnEIC_norm_fun <- function(x, y) {
    return(sqrt(sum(unlist(x) * unlist(y))))
  }

  PnEIC_wt_fun <- function(PnEIC) {
    return(PnEIC)
  }


  # initial estimation + SL / HAL ----------------------------------------------------------

  ## cross validation setup ----
  # stratifying cv so that folds are balanced for treatment assignment & outcomes
  StrataIDs <- as.numeric(factor(paste0(Data$Trt, ":", Data$Event)))
  CVFolds <- origami::make_folds(Data, origami::folds_vfold, strata_ids = StrataIDs)

  ## PropScore score ----
  TrtTask <- sl3::make_sl3_Task(
    data = Data[, -c("Time", "Event", "ID")],
    covariates = colnames(Data[, -c("Time", "Event", "ID", "Trt")]),
    outcome = "Trt"
  )

  TrtSL <- sl3::Lrnr_sl$new(learners = Models[["A"]],
                            folds = CVFolds)
  TrtFit <- TrtSL$train(TrtTask)
  PropScore <- TrtFit$predict()

  if (min(PropScore) < PropScoreCutoff | max(PropScore) > 1 - PropScoreCutoff) {
    warning("practical positivity violation likely, truncating PropScore +-0.05\n")
    PropScore[PropScore < 0.05] <- 0.05
    PropScore[PropScore > 0.95] <- 0.95
  }

  ## hazards: Events & censoring ----
  SupLrnModels <- list()
  for (j in grep("\\d+", names(Models), value = T)) {
    Models_j <- Models[[j]]
    SupLrnLibRisk <- data.table(matrix(NaN, nrow = nrow(Data), ncol = length(Models_j)))
    colnames(SupLrnLibRisk) <- names(Models_j)

    for (Fold_v in CVFolds) {
      TrainIndices <- Fold_v$training_set
      ValidIndices <- Fold_v$validation_set
      TrainData <- Data[TrainIndices, -c("ID")]
      ValidData <- Data[ValidIndices, ][order(-Time)]

      for (k in 1:length(Models_j)) {
        ## train model ----
        CoxphArgs <- list("formula" = Models_j[[k]], "data" = TrainData)
        ModelFit <- do.call(coxph, CoxphArgs)

        ## validation loss (-log partial likelihood) ----
        ValidData[, fit.lp := predict(ModelFit, type = "lp", newdata = ValidData)]
        ValidData[, at.risk := cumsum(exp(fit.lp))]
        ValidData[at.risk == 0, at.risk := 1]
        ValidData[, names(Models_j)[k] := (Event == j) * (fit.lp - log(at.risk))]
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
    ModelFit <- do.call(coxph, list("formula" = SLMod$SupLrnModel, "data" = Data))
    Hazards <- rbind(data.table(time = 0, hazard = 0),
                     suppressWarnings(setDT(basehaz(ModelFit, centered = TRUE))))
    colnames(Hazards) <- c("Time", "SupLrnHaz")
    return(list("fit" = ModelFit, "bhaz" = Hazards))
  })

  ## output initial estimates ----
  HazTimes <- unique(c(TargetTimes, Data$Time))
  HazTimes <- HazTimes[HazTimes <= max(TargetTimes)]

  Hazards <- data.table("Time" = c(0, HazTimes))[order(Time)]
  Hazards <- rbind(cbind(Hazards, A = 0), cbind(Hazards, A = 1))

  ## baseline hazards ----
  for (j in grep("\\d+", names(Models), value = T)) {
    Hazards <- merge(Hazards, SupLrnFits[[j]][["bhaz"]], by = "Time", all.x = T)
    Hazards[, SupLrnHaz := zoo::na.locf(SupLrnHaz), by = "A"]
    if (j == "0") {
      Hazards[, SupLrnHaz := c(0, SupLrnHaz[-.N]), by = "A"]
      setnames(Hazards, "SupLrnHaz", "cumhaz.C.t-")
    } else {
      # hazards[, paste0("cumbhaz.j", j) := haz]
      Hazards[, SupLrnHaz := c(0, diff(SupLrnHaz)), by = "A"]
      setnames(Hazards, "SupLrnHaz", paste0("bhaz.j", j))
    }
  }

  ## clever covariate data.table Ht
  Ht <- rbind(cbind(Data, A = 1),
              cbind(Data, A = 0))
  setcolorder(Ht, c("ID", "Trt", "Time", "Event"))

  ## clever covariate A component ----
  Ht[A == 1, Ht.a := 1 / PropScore]
  Ht[A == 0, Ht.a := 1 / (1 - PropScore)]

  ## clever covariate hazard components ----

  setnames(Ht, c("A", "Trt"), c("Trt", "A"))
  for (j in grep("\\d+", names(Models), value = T)) {
    if (j == "0") {
      Ht[, cox.C := predict(SupLrnFits[["0"]]$fit, newdata = Ht, type = "risk")]
    } else {
      Ht[, (paste0("cox.j", j)) := predict(SupLrnFits[[j]]$fit,
                                           newdata = Ht, type = "risk")]
    }
  }
  setnames(Ht, c("A", "Trt"), c("Trt", "A"))
  # A = target trt, Trt = assigned trt; Time = eval time, T_tilde = obs event time
  setnames(Ht, "Time", "T_tilde")

  Ht_colnames <- c("ID", "Trt", "A", "T_tilde", "Event", "Ht.a",
                   grep("cox", colnames(Ht), value = T))
  Ht <- Ht[, ..Ht_colnames]

  Ht <- merge(Ht, Hazards, by = "A", allow.cartesian = T)[order(ID, Time, A)]
  setcolorder(Ht, c("ID", "Time", "A"))

  ### clever covariate C component ----
  Ht[, "S.C.a.t-" := exp(-`cumhaz.C.t-` * cox.C)]
  if (min(Ht$`S.C.a.t-`) < 0.05) {
    warning(paste0("targeting a time where probability of being censored >95%,",
                   "cens survival truncated at min 0.05, but you should really",
                   "aim lower\n"))
  }
  Ht[, "Ht.g" := Ht.a / `S.C.a.t-`]
  Ht <- Ht[, -c("Ht.a", "cox.C", "cumhaz.C.t-", "S.C.a.t-")]

  ### event-free survival ----
  Ht[, S.t := 1]
  ## expand out bhaz and cox to js, then do by = c("ID", "A", "J")
  for (j in TargetEvents) {
    Ht[, S.t := S.t * exp(-cumsum(get(paste0("bhaz.j", j)) *
                                    get(paste0("cox.j", j)))),
       by = c("ID", "A")]
  }
  if (min(Ht[, S.t]) <= 0) {
    stop("you're targeting a time when there are no survivors, aim lower\n")
  }
  Ht[, "S.t-" := c(1, S.t[-.N]), by = c("ID", "A")]

  ### cause-specific risks ----
  for (j in TargetEvents) {
    Ht[, (paste0("F.j", j, ".t")) := cumsum(`S.t-` *
                                              get(paste0("bhaz.j", j)) *
                                              get(paste0("cox.j", j))
    ),
    by = c("ID", "A")]
  }


  ### cause-specific risks at target times ----
  ### target event (j), target time (k), cause (l) clever covariate ----
  for (k in TargetTimes) {
    for (j in TargetEvents) {
      #### event j risk at target time k ----
      Ht[, fjtau := get(paste0("F.j", j, ".t"))[Time == k], by = c("ID", "A")]
      for (l in Events) {
        Ht[, hjlt := 1] # why was 1 the default, a 0 should squash oversights no?
        #### event component ----
        Ht[, hjlt := (l == j) - (fjtau - get(paste0("F.j", j, ".t"))) / S.t]
        setnames(Ht, "hjlt", paste0("h.j", j, ".l", l, ".t", k))
      } #### mean treated/control event j risk at target time k ----
      Ht[, paste0("Psi.j", j, ".t", k) := mean(fjtau), by = "A"]
      setnames(Ht, "fjtau", paste0("F.j", j, ".t", k))
    }
  }

  ## initial estimator (g-computation) ----
  Psi_init <- Ht[, mget(c("A", grep("Psi.+", colnames(Ht), value = T)))]
  Psi_init <- Psi_init[, lapply(.SD, unique), by = "A"]
  Psi_init <- melt(Psi_init, id.vars = "A")
  Psi_init[, c("dummy", "event", "time") := tstrsplit(variable, "\\.(j|t)")]
  Psi_init <- Psi_init[, -c("variable", "dummy")]
  Psi_init <- dcast(Psi_init, A + time ~ event, value.var = "value")
  Psi_init[, S := do.call(sum, mget(paste0(1:3))), by = c("A", "time")]
  Psi_init <- Psi_init[, lapply(.SD, as.numeric)][order(time, A)]


  # Update step ---------------------------------------------------------------------------------

  ## EIC ----
  Ht_eic <- data.table(ID = Data$ID)
  for (k in TargetTimes) {
    for (j in TargetEvents) {
      Ht[, EIC := 0]
      for (l in Events) {
        Ht[, EIC := EIC +
             Ht.g * (Time <= k) * (Time <= T_tilde) * (Trt == A) *
             get(paste0("h.j", j, ".l", l, ".t", k)) *
             ((Time == T_tilde ) * (Event == l) -
                get(paste0("cox.j", l)) * get(paste0("bhaz.j", l)))]
      }
      Ht_eic <- merge(Ht_eic,
                      dcast(Ht[, .(EIC = sum(EIC),
                                   Fjt = unique(get(paste0("F.j", j, ".t", k))),
                                   Psijt = unique(get(paste0("Psi.j", j, ".t", k)))),
                               by = c("ID", "A")][, EIC := EIC + Fjt - Psijt],
                            ID ~ A, value.var = "EIC"),
                      by = "ID")
      setnames(Ht_eic, "0", paste0("EIC.a0.j", j, ".t", k))
      setnames(Ht_eic, "1", paste0("EIC.a1.j", j, ".t", k))
    }
  }

  init_ic <- list("ic" = Ht_eic[, -c("ID")],
                  "pnEIC" = colMeans(Ht_eic[, -c("ID")]),
                  "seEIC" = sqrt(diag(var(Ht_eic[, -c("ID")]))))

  PnEIC <- init_ic$pnEIC
  PnEIC <- as.data.table(cbind("A" = 1:0,
                               rbind(PnEIC[grep("a1", names(PnEIC))],
                                     PnEIC[grep("a0", names(PnEIC))])))
  setnames(PnEIC, 2:ncol(PnEIC),
           do.call(paste0, expand.grid("PnEIC.j", TargetEvents, ".t", TargetTimes)))

  PnEIC_wtd <- PnEIC_wt_fun(PnEIC)

  Ht <- merge(Ht, PnEIC_wtd, by = "A")
  PnEIC_norm <- PnEIC_norm_fun(PnEIC, PnEIC_wtd)

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
    for (j in TargetEvents) { # loop over target events
      for (k in TargetTimes) { # loop over target times
        for (l in Events) { # nested loop over all events
          Ht[, (paste0("h.j", j, ".l", l, ".t", k, ".tmp")) :=
               get(paste0("h.j", j, ".l", l, ".t", k))]
        }
      }
      Ht[, (paste0("cox.j", j, ".tmp")) := get(paste0("cox.j", j))]
    }

    ## 5.1 calculate update step direction -------------------------------------

    for (l in Events) { # loop over event hazards
      Ht[, (paste0("delta.j", l, ".dx")) := 0] # fluctuation of causes-specific hazard
      for (j in TargetEvents) { # loop over target Events
        for (k in TargetTimes) { # loop over target times
          Ht[, (paste0("delta.j", l, ".dx")) := get(paste0("delta.j", l, ".dx")) +
               (Time <= k) * Ht.g *
               get(paste0("h.j", j, ".l", l, ".t", k)) *
               get(paste0("PnEIC.j", j, ".t", k)) / PnEIC_norm]
        }
      }
    }

    ## 5.3 update cause specific hazards ----------------------------------------
    Ht[, S.t := 1]
    for (l in Events) { # loop over competing Events
      Ht[, (paste0("cox.j", l)) := get(paste0("cox.j", l)) *
           exp(OneStepEps * get(paste0("delta.j", l, ".dx")))]
      #### truncation cox fit +- 500 why ? -----
      Ht[, (paste0("cox.j", l)) := pmax(-500, pmin(500, get(paste0("cox.j", l))))]
      # Does changing the order of event updates matter? -------------------------------
      Ht[, S.t := S.t * exp(-cumsum(get(paste0("bhaz.j", l)) *
                                      get(paste0("cox.j", l)))),
         by = c("ID", "A")]
      # when does cox go to Infinity? -------
      if (any(is.infinite(Ht[, get(paste0("cox.j", l))])))
        stop(paste0("cox.fit for j=", j, " went to infinity. wtf happened?\n"))
      if (min(Ht[, S.t]) <= 0)
        stop(paste0("Survival at some time is leq 0. wtf happened?\n"))
    }
    Ht[, "S.t-" := c(1, S.t[-.N]), by = c("ID", "A")]

    ## update clever covariates ----
    for (j in TargetEvents) {
      Ht[, (paste0("F.j", j, ".t")) := cumsum(`S.t-` *
                                                get(paste0("bhaz.j", j)) *
                                                get(paste0("cox.j", j))),
         by = c("ID", "A")]
      for (k in TargetTimes) {
        # Ht[, (paste0("S.t", k)) := S.t[Time == k], by = c("ID", "A")]
        Ht[, fjtau := get(paste0("F.j", j, ".t"))[Time == k], by = c("ID", "A")]
        for (l in Events) {
          Ht[, hjlt := 1] #### ASK HELENE: why is 1 the default? ---------
          #### event component ----
          Ht[, hjlt := (l == j) - (fjtau - get(paste0("F.j", j, ".t"))) / S.t]

          Ht[round(S.t, 8) == 0, hjlt := 0 - (l == j)] ##### ASK HELENE WHY -------

          Ht[, (paste0("h.j", j, ".l", l, ".t", k)) := hjlt]
        } #### mean treated/control event j risk at target time k ----
        Ht[, (paste0("Psi.j", j, ".t", k)) := mean(fjtau), by = "A"]
        Ht[, (paste0("F.j", j, ".t", k)) := fjtau]
      }
    }
    Ht <- Ht[, -c("hjlt", "fjtau")]
    ## recalculate eic -----
    Ht_eic <- data.table(ID = Data$ID)
    for (k in TargetTimes) {
      for (j in TargetEvents) {
        Ht[, EIC := 0]
        for (l in Events) {
          Ht[, EIC := EIC + get(paste0("h.j", j, ".l", l, ".t", k)) * Ht.g *
               (Time <= k) * (Time <= T_tilde) * (Trt == A) *
               ((Time == T_tilde) * (Event == l) -
                  get(paste0("cox.j", l)) * get(paste0("bhaz.j", l)))]
        }
        Ht_eic <- merge(Ht_eic,
                        dcast(Ht[, .(EIC = sum(EIC),
                                     Fjt = unique(get(paste0("F.j", j, ".t", k))),
                                     Psijt = unique(get(paste0("Psi.j", j, ".t", k)))),
                                 by = c("ID", "A")][, EIC := EIC + Fjt - Psijt],
                              ID ~ A, value.var = "EIC"),
                        by = "ID")
        setnames(Ht_eic, "0", paste0("EIC.a0.j", j, ".t", k))
        setnames(Ht_eic, "1", paste0("EIC.a1.j", j, ".t", k))
      }
    }

    ic <- list("ic" = Ht_eic,
               "pnEIC" = colMeans(Ht_eic[, -c("ID")]),
               "seEIC" = sqrt(diag(var(Ht_eic[, -c("ID")]))))

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
      for (j in TargetEvents) { # loop over target Events
        for (k in TargetTimes) { # loop over target times
          for (l in Events) { # nested loop over Events
            Ht[, (paste0("h.j", j, ".l", l, ".t", k)) :=
                 get(paste0("h.j", j, ".l", l, ".t", k, ".tmp"))]
          }
        }
        Ht[, (paste0("cox.j", j)) := get(paste0("cox.j", j, ".tmp"))]
      }

      PnEIC <- PnEIC_prev
      PnEIC_wtd <- PnEIC_wtd_prev
      PnEIC_norm <- PnEIC_norm_prev
    } else {
      # update PnEIC columns in Ht with updated values
      ht_eic_cols <- grep("PnEIC", colnames(Ht), value = T)
      Ht <- Ht[, !..ht_eic_cols]
      Ht[PnEIC, on = .(A = A), (ht_eic_cols) := mget(sprintf("i.%s", ht_eic_cols))]
      gc() # otherwise memory blows up

      summ_eic <- melt(ic$ic, id.vars = "ID")
      summ_eic[, c("dummy", "A", "J", "T") := tstrsplit(variable, "\\.(a|j|t)")]
      summ_eic <- summ_eic[, -c("variable", "dummy")]
      summ_eic[, A := paste0("a", get("A"))]
      summ_eic[, J := paste0("j", get("J"))]
      summ_eic[, "T" := paste0("t", get("T"))]
      summ_eic <- dcast(summ_eic, ID + A + `T` ~ J, value.var = "value")[
        , S := do.call(sum, mget(paste0("j", 1:3))),
        by = c("ID", "A", "T")]
      summ_eic <- dcast(summ_eic, ID ~ ...,
                        value.var = c("j1", "j2", "j3", "S"),
                        sep = ".")

      onestep_stop <- abs(colMeans(summ_eic[, -"ID"])) <=
        sqrt(colMeans(summ_eic[, -"ID"]^2)) / ( sqrt(nrow(Data)) * log(nrow(Data)))

      if (Verbose)
        cat(signif(abs(colMeans(summ_eic[, -"ID"])) /
                     (sqrt(colMeans(summ_eic[, -"ID"]^2)) /
                     ( sqrt(nrow(Data)) * log(nrow(Data)))), 2), "\n")


      if (all(onestep_stop) | step == NumUpdateSteps) {
        if (step == NumUpdateSteps) {
          message("Warning: Algorithm did not converge")
        }

        if (Verbose) {
          print(paste0("converged", " at ", step, "th step"))
          print(paste0("PnEIC = ", colMeans(summ_eic[, -"ID"])))
        }

        ## format output ----
        Psi_tmle <- Ht[, mget(c("A", grep("Psi.+", colnames(Ht), value = T)))]
        Psi_tmle <- Psi_tmle[, lapply(.SD, unique), by = "A"]
        Psi_tmle <- melt(Psi_tmle, id.vars = "A")
        Psi_tmle[, c("dummy", "event", "time") := tstrsplit(variable, "\\.(j|t)")]
        Psi_tmle <- Psi_tmle[, -c("variable", "dummy")]
        Psi_tmle <- dcast(Psi_tmle, A + time ~ event, value.var = "value")
        Psi_tmle[, S := do.call(sum, mget(paste0(1:3))), by = c("A", "time")]
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



