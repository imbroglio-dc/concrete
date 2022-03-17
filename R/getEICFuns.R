# get EICs
getEIC <- function(PropScores, SupLrnFits, Hazards, TargetEvents, TargetTimes, Events,
                   RegsOfInterest, MinNuisanceDenom, GComp = F) {
    Targets <- expand.grid("Time" = TargetTimes, "Event" = TargetEvents)

    ICs <- lapply(names(RegsOfInterest), function(RegName) {
        PredFits <- getPredFits(RegName, RegsOfInterest, PropScores, Data, SupLrnFits)

        Censored <- !is.null(Hazards[["LagCumHaz.C.t"]])
        # loop over individuals
        IC.a <- lapply(1:nrow(PredFits), function(i) {
            Subj <- PredFits[i, ]
            SubjRisk <- getSubjRisk(Subj, Hazards, MinNuisanceDenom, Events, Targets, Censored)
            SubjIC.a <- getSubjICs(Subj, SubjRisk, Hazards, MinNuisanceDenom, Events, Targets, Censored)
            # browser()
            # if (GComp) {
            #     GCompEst <- getGComp(SubjRisk, TargetTimes)
            #     return(list("IC" = SubjIC.a, "GComp" = GCompEst))
            # } else {
            #     return(list("IC" = SubjIC.a))
            # }
        })

        output <- list()
        IC.a <- as.data.table(do.call(rbind, lapply(1:length(IC.a),
                                                    function(i) cbind("ID" = i, IC.a[[i]]))))
        IC.a[, IC := IC + `F` - mean(`F`), by = c("Time", "Event")]
        tmp <- IC.a[, list(-1, -sum(IC), 1 - sum(`F`)), by = c("Time", "ID")]
        setnames(tmp, c(paste0("V", 1:3)), c("Event", "IC", "F"))
        IC.a <- rbind(IC.a, tmp)

        output[["EIC"]] <- cbind("A" = RegName, IC.a)

        output[["SummEIC"]] <- output[["EIC"]][, list(mean(IC), mean(`F`), sqrt(var(IC)/.N)),
                                               by = c("Time", "Event")]
        setnames(output[["SummEIC"]], paste0("V", 1:3), c("PnEIC", "Psi", "sePsi.IC"))
        output[["SummEIC"]] <- cbind("A" = RegName, output[["SummEIC"]])
        return(output)
    })

    return(list("EIC" = do.call(rbind, lapply(ICs, function(a) a[["EIC"]])),
                "SummEIC" = do.call(rbind, lapply(ICs, function(a) a[["SummEIC"]]))))
}

getPredFits <- function(RegName, RegsOfInterest, PropScores, Data, SupLrnFits) {
    PredFits <- as.data.table(cbind(Data, "A" = RegsOfInterest[[RegName]]))
    PredFits[, g.star := attr(RegsOfInterest[[RegName]], "g.star")]

    ## clever covariate hazard components ----
    setnames(PredFits, c("A", "Trt"), c("Trt", "A"))
    for (j in grep("\\d+", names(SupLrnFits), value = T)) {
        if (j == "0") {
            PredFits[, Cox.C := predict(SupLrnFits[["0"]]$fit, newdata = PredFits, type = "risk")]
        } else {
            PredFits[, (paste0("Cox.j", j)) := predict(SupLrnFits[[j]][["fit"]],
                                                       newdata = PredFits, type = "risk")]
        }
    }
    # A = target trt, Trt = assigned trt; Time = eval time, T.tilde = observed event time
    setnames(PredFits, c("A", "Trt", "Time"), c("Trt", "A", "T.tilde"))

    PredFitsColNames <- c("ID", "Trt", "A", "g.star","T.tilde", "Event",
                          grep("Cox", colnames(PredFits), value = T))
    PredFits <- PredFits[, PredFitsColNames, with = FALSE]

    PredFits[, Ht.a := 1 / PropScores[[RegName]]]
    return(PredFits[])
}

getSubjRisk <- function(Subj, Hazards, MinNuisanceDenom, Events, Targets, Censored) {
    ### clever covariate nuisance component ----
    Ht.g <- Subj[["Ht.a"]]
    if (Censored)
        Ht.g <- Ht.g / exp(-Hazards[["LagCumHaz.C.t"]] * Subj[["Cox.C"]])

    truncated <- FALSE
    if (max(Ht.g) > 1 / MinNuisanceDenom) {
        Ht.g <- 1 / truncNuisanceDenom(1/Ht.g, MinNuisanceDenom)
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

getSubjICs <- function(Subj, SubjRisk, Hazards, MinNuisanceDenom, Events, Targets, Censored) {
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

