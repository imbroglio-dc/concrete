# get EICs
getEIC <- function(Estimates, Data, RegsOfInterest, Censored, TargetEvents, TargetTimes, Events, MinNuisance, GComp = F) {
    Targets <- expand.grid("Time" = TargetTimes, "Event" = TargetEvents)
    EvalTimes <- attr(Estimates, "times")
    T.tilde <- Data[["Time"]]
    Delta <- Data[["Event"]]

    ICs <- lapply(Estimates, function(RegimeEsts) {
        NuisanceWeight <- sapply(1:length(RegimeEsts[["PropScore"]]), function(i) {
            RegimeEsts[["PropScore"]][i] * RegimeEsts[["Survival"]][["LaggedCensSurv"]][, i]})
        NuisanceWeight <- 1 / truncNuisanceDenom(NuisanceWeight, MinNuisance)
        GStar <- attr(PropScore, "g.star")
        Hazards <- RegimeEsts[["Hazards"]]
        TotalSurv <- RegimeEsts[["Survival"]][["TotalSurv"]]

        IC.a <- getICs(GStar, Hazards, TotalSurv, NuisanceWeight, Targets, Events, EvalTimes, GComp)

        getICs <- function(GStar, Hazards, TotalSurv, NuisanceWeight, Targets,
                           Events, EvalTimes, GComp) {

            IC.a <- lapply(1:length(PropScore), function(i) {
                Hazards.i <- lapply(Hazards, function(haz) haz[, i])
                Surv.i <- TotalSurv[, i]
                Risks.i <- lapply(Hazards, function(haz) cumsum(Surv.i * haz))
                LagSc.i <- LaggedCensSurv[, i]

                if (GStar[i] == 0) # 1(A == a*)
                    return(cbind(Targets, "IC" = 0,
                                 "F" = apply(Targets, 1, function(target) {
                                     tau <- target[["Time"]]
                                     j <- target[["Event"]]
                                     return(Risks.i[[as.character(j)]][EvalTimes == tau])
                                 })))

                IC.jk <- t(apply(Targets, 1, function(target) {
                    t.tilde <- T.tilde[i]
                    tau <- target[["Time"]]
                    time.cutoff <- min(tau, t.tilde) ## 1(t \leq tau) * 1(t \leq t.tilde)
                    j <- target[["Event"]]
                    F.j.tau <- Risks.i[[as.character(j)]][EvalTimes == tau]

                    Surv.ik <- Surv.i[EvalTimes <= time.cutoff]
                    haz.J.ik <- lapply(Hazards.i, function(r) r[EvalTimes <= time.cutoff])
                    F.j.t <- Risks.i[[as.character(j)]][EvalTimes <= time.cutoff]

                    IC.jk <- sum(sapply(Events, function(l) {
                        h.jk <- GStar[i] * NuisanceWeight[i] * ((l == j) - (F.j.tau - F.j.t) / Surv.ik)

                        IC.ljk <- sum(h.jk * ((EvalTimes == t.tilde) * (Delta[i] == l) - haz.J.ik[[as.character(l)]]))
                        ## the second component ( ... + F_j(tau | a, L) - Psi ) is done outside
                        return(IC.ljk)
                    }))
                    return(c("IC" = IC.jk, "F" = F.j.tau))
                }))
                IC.jk <- cbind(Targets, IC.jk)
                return(IC.jk)
            })
            return(list("IC" = IC.a, "NuisanceWeights" = NuisanceWeights))
        }

        # loop over individuals
        IC.a <- lapply(1:length(PropScore), function(i) {
            SubjPropScore <- PropScore[i]
            Subjg.star <- g.star[i]
            SubjHazards <- lapply(Hazards, function(haz) haz[, i])
            SubjTotalSurv <- TotalSurv[, i]
            SubjLaggedCensSurv <- LaggedCensSurv[, i]

            SubjRisk <- getSubjRisk(SubjPropScore, SubjHazards, SubjTotalSurv, SubjLaggedCensSurv,
                                    MinNuisance, Events, Targets, Censored)

            SubjIC.a <- getSubjICs(SubjPropScore, SubjHazards, SubjTotalSurv, SubjLaggedCensSurv,
                                   MinNuisance, Events, Targets, Censored)
            return(list("IC" = SubjIC.a, "Ht.g" = SubjRisk[["Ht.g"]], "S.t" = SubjRisk[["S.t"]]))
        })

        output <- list()
        output[["PredFits"]] <- PredFits
        output[["Ht.g"]] <- as.data.table(do.call(cbind, lapply(IC.a, function(i) i[["Ht.g"]])))
        output[["S.t"]] <- as.data.table(do.call(cbind, lapply(IC.a, function(i) i[["S.t"]])))

        IC.a <- as.data.table(do.call(rbind, lapply(1:length(IC.a),
                                                    function(i) cbind("ID" = i, IC.a[[i]][["IC"]]))))
        IC.a[, IC := IC + `F` - mean(`F`), by = c("Time", "Event")]
        tmp <- IC.a[, list(-1, -sum(IC), 1 - sum(`F`)), by = c("Time", "ID")]
        setnames(tmp, c(paste0("V", 1:3)), c("Event", "IC", "F"))
        output[["EIC"]] <- rbind(IC.a, tmp)

        output[["SummEIC"]] <- output[["EIC"]][, list(mean(IC), mean(`F`), sqrt(var(IC)/.N)),
                                               by = c("Time", "Event")]
        setnames(output[["SummEIC"]], paste0("V", 1:3), c("PnEIC", "Psi", "sePsi.IC"))

        return(output)
    })
            names(ICs) <- names(RegsOfInterest)

            out <- list("EIC" = lapply(ICs, function(a) a[["EIC"]]),
                        "SummEIC" = lapply(ICs, function(a) a[["SummEIC"]]),
                        "PredFits" = lapply(ICs, function(a) a[["PredFits"]]))
            return(ICs)
}


getSubjRisk <- function(Subj, Hazards, MinNuisance, Events, Targets, Censored) {
    ### clever covariate nuisance component ----
    Ht.g <- Subj[["Ht.a"]]
    if (Censored)
        Ht.g <- Ht.g / exp(-Hazards[["LagCumHaz.C.t"]] * Subj[["Cox.C"]])

    truncated <- FALSE
    if (max(Ht.g) > 1 / MinNuisance) {
        Ht.g <- 1 / truncNuisanceDenom(1/Ht.g, MinNuisance)
        truncated <- TRUE
    }

    ### event-free survival ----
    S.t <- exp(-cumsum(rowSums(sapply(Events, function(j) {
        Hazards[[paste0("BaseHaz.j", j)]] * Subj[[paste0("Cox.j", j)]]
    }))))
    # LagS.t <- c(1, head(S.t, -1))

    ### cause-specific risks ----
    SubjRisk <- do.call(cbind, lapply(Events, function(j) {
        cumsum(S.t * Hazards[[paste0("BaseHaz.j", j)]] * Subj[[paste0("Cox.j", j)]])
    }))
    SubjRisk <- as.data.table(do.call(cbind, list(Hazards[["Time"]], Ht.g, S.t, SubjRisk)))
    colnames(SubjRisk) <- c("Time", "Ht.g", "S.t", paste0("F.j", Events, ".t"))
    ## SURVIVAL + RISKS DON'T SUM TO 1 - because naive implementation, fix? ----
    return(SubjRisk)
}

getGComp <- function(SubjRisk, Targets) {
    RiskCols <- c("Time", grep("F", names(SubjRisk), value = TRUE))
    GCompEst <- SubjRisk[Time %in% Targets[["Time"]], .SD, .SDcols = RiskCols]
    GCompEst <- melt()
    return(SubjRisk[Time %in% TargetTimes, -c("Ht.g")])
}

getSubjICs <- function(SubjPropScore, SubjHazards, SubjTotalSurv, SubjLaggedCensSurv,
                       MinNuisance, Events, Targets, Censored) {
    if (Subj[["g.star"]] == 0) # 1(A == a*)
        return(cbind(Targets, "IC" = 0,
                     "F" = apply(Targets, 1, function(target) {
                         tau <- target[["Time"]]
                         j <- target[["Event"]]
                         return(SubjRisk[[paste0("F.j", j, ".t")]][SubjRisk[["Time"]] == tau])
                     })))

    IC.jk <- t(apply(Targets, 1, function(target) {
        t.tilde <- Subj[["T.tilde"]]
        tau <- target[["Time"]]
        j <- target[["Event"]]
        F.j.tau <- SubjRisk[[paste0("F.j", j, ".t")]][SubjRisk[["Time"]] == tau]
        SubjRisk.ik <- SubjRisk[Time <= min(tau, t.tilde), ] ## 1(t \leq tau) * 1(t \leq t.tilde)
        Hazards.ik <- Hazards[Time <= min(tau, t.tilde), ] ## 1(t \leq tau) * 1(t \leq t.tilde)
        F.j.t <- SubjRisk.ik[[paste0("F.j", j, ".t")]]

        IC.jk <- sum(sapply(Events, function(l) {
            h.jk <- Subj[["g.star"]] * SubjRisk.ik[["Ht.g"]] *
                ((l == j) - (F.j.tau - F.j.t) / SubjRisk.ik[["S.t"]])

            IC.ljk <- sum(h.jk * ((SubjRisk.ik[["Time"]] == t.tilde) * (Subj[["Event"]] == l) -
                                      Hazards.ik[[paste0("BaseHaz.j", l)]] * Subj[[paste0("Cox.j", l)]]))
            ## the second component ( ... + F_j(tau | a, L) - Psi ) is done outside
            return(IC.ljk)
        }))
        return(c("IC" = IC.jk, "F" = F.j.tau))
    }))
    IC.jk <- cbind(Targets, IC.jk)
    return(IC.jk)
}

