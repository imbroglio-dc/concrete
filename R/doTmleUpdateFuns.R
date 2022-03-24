# Pseudocode
# Variables:
#     N = number of observations
#     T = number of time points
#     J = number of target Events
#     L = number of Events
#     K = number of target time points
# Data:
#     For a in A:
#         pi(A | W)     : N x 1 vector
#         Sc(t | A, W)  : T x N matrix
# updated S(t | A, W)   : T x N matrix
# updated h(t | A, W)   : L (T x N) matrices
# Do:
#     For a in A:
#         for l in L:
#             for j in J and tau in K:
#                 for i in N:
#                     calculate       h_jk      = (T x (J x K)) matrix          # clever covariate
#                     fetch           PnEIC     = (J x K) x 1 vector            # emprcl mean EIC
#                     calculate  <PnEIC, h_jk>  = T x 1 vector                  # update direction
#           store     recalculate  h(t | A, W)  = T x 1 vector                  # updated hazard
#
#           store     recalculate S(t | A, W)   = (T x N) matrices
#
#         for j in J , tau in K: recalculate EIC :
#
#     undo or commit updated h() and F() based on PnEIC

doTmleUpdate <- function(Estimates, SummEIC, Data, Censored, TargetEvents, TargetTimes, Events,
                         NumUpdateSteps, OneStepEps, Verbose) {
    Targets <- expand.grid("Time" = TargetTimes, "Event" = TargetEvents)
    EvalTimes <- attr(Estimates, "times")
    T.tilde <- Data[["Time"]]
    Delta <- Data[["Event"]]

    for (step in 1:NumUpdateSteps) {
        ## Check if EIC is solved sufficienty and return outputs ----
        ## check PnEIC <= seEIC / (sqrt(n) log(n))
        OnestepStop <- SummEIC[, list("check" = abs(PnEIC) <= `seEIC/(root(n)log(n))`,
                                      "ratio" = abs(PnEIC) / `seEIC/(root(n)log(n))`),
                               by = c("Trt", "Time", "Event")]
        if (Verbose)
            print(OnestepStop[["ratio"]])
        if (all(sapply(OnestepStop[["check"]], isTRUE))) {

            return()
        }

        ## one-step tmle loop starts here ----

        if (Verbose)
            cat("starting step", step, "with update epsilon =", OneStepEps, "\n")

        ## make backups ----
        ## save current cause-specific hazards and clever covs in case the update
        ## makes things worse

        PrevSummEIC <- SummEIC
        PrevEsts <- Estimates

        ## Get updated hazards
        updatedHazards <- lapply(Estimates, function(est.a) {
            NuisanceWeight <- RegimeEsts[["NuisanceWeight"]]
            GStar <- attr(RegimeEsts[["PropScore"]], "g.star.intervention")
            Hazards <- RegimeEsts[["Hazards"]]
            TotalSurv <- RegimeEsts[["Survival"]][["TotalSurv"]]

            lapply(Hazards, function() {})
        })
        updateHazards <- function(GStar, Hazards, TotalSurv, NuisanceWeight, Targets,
                           Events, EvalTimes, GComp) {
            # loop over individuals
            lapply(Estimates, function(est.a) {
                tmp <- lapply(est.a[["Hazards"]], function(haz.al) {
                    haz.ali <- sapply(1:ncol(est.a[["NuisanceWeight"]]), function(i) {

                    })
                })
            })

            IC.a <- do.call(rbind, lapply(1:length(PropScore), function(i) {
                LagSc.i <- LaggedCensSurv[, i]
                Surv.i <- TotalSurv[, i]
                Hazards.i <- lapply(Hazards, function(haz) haz[, i])
                Risks.i <- lapply(Hazards.i, function(haz.i) cumsum(Surv.i * haz.i))

                if (GStar[i] == 0) # 1(A == a*)
                    return(cbind("ID" = i, Targets, "IC" = 0,
                                 "F.j.tau" = apply(Targets, 1, function(target) {
                                     tau <- target[["Time"]]
                                     j <- target[["Event"]]
                                     return(Risks.i[[as.character(j)]][EvalTimes == tau])
                                 })))

                IC.jk <- t(apply(Targets, 1, function(target) {
                    j <- target[["Event"]]
                    tau <- target[["Time"]]
                    t.tilde <- T.tilde[i]
                    time.cutoff <- min(tau, t.tilde) ## 1(t \leq tau) * 1(t \leq t.tilde)
                    F.j.tau <- Risks.i[[as.character(j)]][EvalTimes == tau]
                    s.ik <- EvalTimes[EvalTimes <= time.cutoff]
                    Surv.ik <- Surv.i[EvalTimes <= time.cutoff]
                    haz.J.ik <- lapply(Hazards.i, function(r) r[EvalTimes <= time.cutoff])
                    F.j.t <- Risks.i[[as.character(j)]][EvalTimes <= time.cutoff]

                    IC.jk <- sum(sapply(Events, function(l) {
                        h.jk <- GStar[i] * NuisanceWeight[i] * ((l == j) - (F.j.tau - F.j.t) / Surv.ik)
                        IC.ljk <- sum(h.jk * ((s.ik == t.tilde) * (Delta[i] == l) - haz.J.ik[[as.character(l)]]))
                        ## the second EIC component ( ... + F_j(tau | a, L) - Psi ) is done outside
                        return(IC.ljk)
                    }))
                    return(c("IC" = IC.jk, "F.j.tau" = F.j.tau))
                }))
                IC.jk <- cbind("ID" = i, Targets, IC.jk)
                return(IC.jk)
            }))

            ## the second EIC component ( ... + F_j(tau | a, L) - Psi )
            IC.a <- as.data.table(IC.a)
            IC.a[, IC := IC + F.j.tau - mean(F.j.tau), by = c("Time", "Event")]
            return(IC.a[,.SD, .SDcols = !"F.j.tau"])
        }

        UpdatedHazards <-

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
}


## get intervention clever covariates for update ----
getClevCovs <- function(Estimates, Data, RegsOfInterest, Censored, TargetEvents, TargetTimes, Events) {
    if (!all(sapply(list(InitEIC, Hazards), is.null))) {

    } else if (!all(sapply(list(PredFits, Hazards), is.null))) {

    } else {
        warning("Argument(s) missing, gotta provide all update args or all initi args")
    }
}



