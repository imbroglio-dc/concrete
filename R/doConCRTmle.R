
#' Title
#'
#' @param EventTime : Numeric vector (N x 1)
#' @param EventType : Numeric (integer) vector (N x 1)
#' @param Treatment : Numeric vector (N x 1)
#' @param Intervention : list of function (length = A*)
#' @param CovDataTable : data.table (N x ?)
#' @param ID : vector (N x 1)
#' @param TargetTimes : numeric vector (length = K)
#' @param TargetEvents : numeric vector \subset EventType (length = J)
#' @param Models : list of functions (length = L)
#' @param CVArgs : list
#' @param NumUpdateSteps : numeric
#' @param OneStepEps : numeric
#' @param MinNuisance : numeric
#' @param Verbose : boolean
#' @param Censored : boolean
#' @param PropScoreBackend : character
#' @param GComp : boolean
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
doConCRTmle <- function(EventTime, EventType, Treatment, Intervention, CovDataTable,
                        ID = NULL, TargetTimes = sort(unique(EventTime)),
                        TargetEvents = NULL, Models, CVArgs = NULL, NumUpdateSteps = 25,
                        OneStepEps = 0.1, MinNuisance = 0.05, PropScoreBackend = "sl3",
                        Verbose = FALSE, GComp = FALSE, ...)
{
  # check & format parameters  ----------------------------------------------------------------
  # args <- checkConCRTmleArgs(EventTime, EventType, Treatment, CovDataTable,
  #                            ID, TargetTimes, TargetEvents, Models, CVArgs, NumUpdateSteps,
  #                            OneStepEps, Verbose)
  Data <- data.table("ID" = ID,
                     "Time" = EventTime,
                     "Event" = EventType,
                     "Trt" = Treatment,
                     CovDataTable)

  Events <- sort(unique(Data$Event))
  Censored <- 0 %in% Events
  Events <- Events[Events > 0]

  # For user input: data + formula ----
  # try use prodlim::EventHistory.frame
  # x=EventHistory.frame(Hist(time,Status,cens.code="censored")~age+sex+intervention(trt)+stage,data=pbc,specials="intervention")
  # names(x)

  # OBS: learners can either work with the original data or with the design matrix (dummies)

  if (is.list(Intervention)) {
    RegsOfInterest <- lapply(Intervention, function(intervene) {
      if (is.function(intervene)) {
        Regime <- do.call(intervene, list(Treatment, CovDataTable))
        if (is.null(attr(Regime, "g.star"))) {
          attr(Regime, "g.star") <- function(a) as.numeric(a == Regime)
          warning("no g.star input, defaulting to the indicator that observed Treatment == desired RegName")
        }
        return(Regime)
      }
      else stop("Intervention must be a list of functions. See doConCRTmle documentation")
    })
  }

  # helper functions  -------------------------------------------------------------------------

  PnEICNorm_fun <- function(x, y) {
    return(sqrt(sum(unlist(x) * unlist(y))))
  }

  PnEIC_wt_fun <- function(PnEIC, Sigma = NULL) {
    ## INCOMPLETE - work needed to match contmle()
    if (!is.null(Sigma)) {
      SigmaInv <- try(solve(Sigma))
      if (any(class(SigmaInv) == "try-error")) {
        SigmaInv <- solve(Sigma + diag(x = 1e-6, nrow = nrow(Sigma)))
        warning("regularization of Sigma needed for inversion")
      }
      PnEIC <- PnEIC %*% SigmaInv
    }
    return(PnEIC)
  }

  # initial estimation ------------------------------------------------------------------------
  InitEsts <- getInitialEstimates(Data, CovdataTable, Models, MinNuisance, TargetEvents,
                                  TargetTimes, RegsOfInterest, PropScoreBackend, Censored)

  # get initial EIC (possibly with GComp Estimate) ---------------------------------------------
  InitEIC <- getEIC(InitEsts, Data, RegsOfInterest, Censored, TargetEvents,
                    TargetTimes, Events, MinNuisance, GComp)
  # EIC <- InitEIC[["EIC"]]
  SummEIC <- InitEIC[["SummEIC"]]
  GCompEst <- InitEIC[["GComp"]]

  ## initial estimator (g-computation) --------------------------------------------------------
  GCompEst <- InitEIC[["SummEIC"]][, c("A", "Time", "Event", "Psi")]

  # Update step -------------------------------------------------------------------------------

  ## EIC ----
  ## one-step tmle loop (one-step) ----

  for (step in 1:NumUpdateSteps) {
    if (Verbose)
      cat("starting step", step, "with update epsilon =", OneStepEps, "\n")
    IC_prev <- IC

    ## make backups ----
    ## save current cause-specific hazards and clever covs in case the update
    ## makes things worse

    OldCols <- grep("(h\\.j)|(Cox\\.j)", colnames(ClevCovTbl), value = TRUE)
    ClevCovTbl.old <- ClevCovTbl[, mget(OldCols)]

    ## 5.1 calculate update step direction -------------------------------------

    getUpdateStepDir <- function(l, TargetEvents, TargetTimes, IC) {
      h.jlk <- expand.grid(paste0("h.j", TargetEvents, ".l", l), paste0(".t", TargetTimes))
      h.jlk <- do.call(paste0, h.jlk)

      ClevCovCols.l <- c("ID", "A", "Time", "Ht.g", h.jlk)
      ClevCov.l <- ClevCovTbl[, .SD, .SDcols = ClevCovCols.l]
      for (a in RegsOfInterest)
        ClevCov.l[A == a, ]

      ClevCov.l <- dcast(ClevCov.l, ... ~ A, value.var = c(h.jlk, "Ht.g"), sep = ".a")
      ClevCov.l <- melt(ClevCov.l, id.vars = c("A", "Time", "Ht.g"))
      ClevCov.l[, c("j", "k") := tstrsplit(variable, "\\.(j|l|t)", keep = c(2, 4))]

      lapply(TargetEvents, function(j) {
        sapply(TargetTimes, function(k) {
          ClevCovCols.ljk <- c("A", "Time", "Ht.g", paste0("h.j", j, ".l", l, ".t", k))
          ClevCovTbl[, list("A" = A, "delta.l.dx" = ), .SDcols = ClevCovCols.ljk]
        })
      })

      for (j in TargetEvents) { # loop over target Events
        for (k in TargetTimes) { # loop over target times
          ClevCovTbl[Time <= k, delta.l := delta.l + Ht.g *
                       get(paste0("h.j", j, ".l", l, ".t", k)) *
                       get(paste0("PnEIC.j", j, ".t", k)) / PnEICNorm]
        }
      }
    }
    ClevCovTbl[, (paste0("delta.j", Events, ".dx ")) := lapply(Events, function(l) 0)]

    for (l in Events) { # loop over event hazards
      ClevCovTbl[, (paste0("delta.j", l, ".dx")) := 0] # fluctuation of causes-specific hazard
      for (j in TargetEvents) { # loop over target Events
        for (k in TargetTimes) { # loop over target times
          ClevCovTbl[Time <= k, (paste0("delta.j", l, ".dx")) := get(paste0("delta.j", l, ".dx")) +
                       Ht.g *
                       get(paste0("h.j", j, ".l", l, ".t", k)) *
                       get(paste0("PnEIC.j", j, ".t", k)) / PnEICNorm]
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

    WtdPnEIC <- PnEIC_wt_fun(PnEIC)
    PnEICNorm <- PnEICNorm_fun(PnEIC, WtdPnEIC)
    if (Verbose)
      cat("Step", step, "mean EIC norm = ", PnEICNorm, "\n")

    if (PnEICNorm_prev <= PnEICNorm) {
      step <- step - 1
      warning(paste0("update overshot! one-step update epsilon of ",
                     OneStepEps, " will be halved\n"))
      OneStepEps <- OneStepEps / 2
      ## Revert to previous cshaz & clevcovs if update made things worse
      ClevCovTbl[, (OldCols) := ClevCovTbl.old[, mget(OldCols)]]

      PnEIC <- PnEIC_prev
      WtdPnEIC <- WtdPnEIC_prev
      PnEICNorm <- PnEICNorm_prev
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
  if (length(Models != length(Events)))
    stop("Models must be provided for every observed event or censoring type")
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

getClevCovs <- function(ClevCovTbl, TargetTimes, TargetEvents, Events) {
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
  SummEIC <- cbind("Summ" = c("PnEIC", "seEIC"),
                   as.data.table(rbind(colMeans(EicTbl[, -c("ID")]),
                                       sqrt(diag(var(EicTbl[, -c("ID")]))))))
  SummEIC <- melt(SummEIC, id.vars = "Summ")
  SummEIC[, c("A", "j", "k") := tstrsplit(variable, "\\.*(a|j|t)", keep = 2:4)]
  SummEIC <- dcast(SummEIC, A + j + k ~ Summ, value.var = "value")

  WtdPnEIC <- PnEIC_wt_fun(SummEIC[["PnEIC"]])
  PnEICNorm <- PnEICNorm_fun(SummEIC[["PnEIC"]], WtdPnEIC)
  SummEIC[, "WtdPnEIC" := WtdPnEIC]

  return(list("IC" = EicTbl[, -c("ID")],
              "SummEIC" = SummEIC,
              "PnEICNorm" = PnEICNorm))
}
